#!/usr/bin/env python3
"""Build country-quiz.html from Natural Earth 1:50m boundary data.

Usage:
    python3 tools/build_quiz.py DATA_DIR [OUT_HTML]

DATA_DIR must contain (from github.com/nvkelso/natural-earth-vector/geojson):
    ne_50m_admin_0_countries.geojson
    ne_50m_admin_0_boundary_lines_land.geojson

For every quiz country this pre-renders one map view: a Lambert azimuthal
equal-area projection centred on that country (so its shape is not stretched),
with neighbouring land, international borders and disputed borders clipped to
the frame. The page then only has to draw ready-made SVG paths, so it works
offline and stays a single file.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point, Polygon, box, shape
from shapely.ops import unary_union

R_KM = 6371.0088
VIEW_W, VIEW_H = 800, 600
MIN_AREA_KM2 = 5000          # smaller countries are left out of the quiz
SIMPLIFY_PX = 0.6

# 193 UN member states (ISO 3166-1 alpha-3), grouped for readability only.
UN_MEMBERS = set("""
DZA AGO BEN BWA BFA BDI CPV CMR CAF TCD COM COD COG CIV DJI EGY GNQ ERI SWZ ETH
GAB GMB GHA GIN GNB KEN LSO LBR LBY MDG MWI MLI MRT MUS MAR MOZ NAM NER NGA RWA
STP SEN SYC SLE SOM ZAF SSD SDN TZA TGO TUN UGA ZMB ZWE
AFG ARM AZE BHR BGD BTN BRN KHM CHN CYP GEO IND IDN IRN IRQ ISR JPN JOR KAZ KWT
KGZ LAO LBN MYS MDV MNG MMR NPL PRK OMN PAK PHL QAT KOR SAU SGP LKA SYR TJK THA
TLS TUR TKM ARE UZB VNM YEM
ALB AND AUT BLR BEL BIH BGR HRV CZE DNK EST FIN FRA DEU GRC HUN ISL IRL ITA LVA
LIE LTU LUX MLT MDA MCO MNE NLD MKD NOR POL PRT ROU RUS SMR SRB SVK SVN ESP SWE
CHE UKR GBR
ATG BHS BRB BLZ CAN CRI CUB DMA DOM SLV GRD GTM HTI HND JAM MEX NIC PAN KNA LCA
VCT TTO USA
ARG BOL BRA CHL COL ECU GUY PRY PER SUR URY VEN
AUS FJI KIR MHL FSM NRU NZL PLW PNG WSM SLB TON TUV VUT
""".split())
assert len(UN_MEMBERS) == 193, len(UN_MEMBERS)
UN_OBSERVERS = {"PSE"}            # Holy See is far below the size cutoff

# Breakaway areas shown as part of the UN-recognised state.
MERGE_INTO = {"SOL": "SOM", "CYN": "CYP"}

# Names as a US middle-school textbook would print them.
NAME_OVERRIDES = {
    "COD": "Democratic Republic of the Congo",
    "COG": "Republic of the Congo",
    "CAF": "Central African Republic",
    "GNQ": "Equatorial Guinea",
    "SSD": "South Sudan",
    "BIH": "Bosnia and Herzegovina",
    "DOM": "Dominican Republic",
    "CIV": "Côte d'Ivoire (Ivory Coast)",
    "SWZ": "Eswatini (Swaziland)",
    "CZE": "Czechia (Czech Republic)",
    "MKD": "North Macedonia",
    "MMR": "Myanmar (Burma)",
    "TLS": "Timor-Leste (East Timor)",
    "USA": "United States",
    "GBR": "United Kingdom",
    "RUS": "Russia",
    "KOR": "South Korea",
    "PRK": "North Korea",
    "LAO": "Laos",
    "VNM": "Vietnam",
    "SYR": "Syria",
    "IRN": "Iran",
    "TZA": "Tanzania",
    "BOL": "Bolivia",
    "VEN": "Venezuela",
    "MDA": "Moldova",
    "BRN": "Brunei",
    "TUR": "Turkey (Türkiye)",
    "BHS": "The Bahamas",
    "GMB": "The Gambia",
    "SLB": "Solomon Islands",
    "ARE": "United Arab Emirates",
    "PSE": "Palestine",
    "CHN": "China",
}

REGION_OF_CONTINENT = {
    "Africa": "africa",
    "Asia": "asia",
    "Europe": "europe",
    "North America": "namerica",
    "South America": "samerica",
    "Oceania": "oceania",
}
MODES = [
    ("africa", "Africa"),
    ("asia", "Asia"),
    ("europe", "Europe"),
    ("namerica", "North & Central America"),
    ("samerica", "South America"),
    ("world", "Whole World"),
]

DASHED_CLASSES = {
    "Disputed (please verify)",
    "Line of control (please verify)",
    "Indefinite (please verify)",
    "Indeterminant frontier",
}


# ---------------------------------------------------------------- projection
def laea(lon, lat, lon0, lat0):
    """Lambert azimuthal equal-area, returns km (x east, y north)."""
    lam = np.radians(lon - lon0)
    phi = np.radians(lat)
    phi0 = math.radians(lat0)
    cosc = math.sin(phi0) * np.sin(phi) + math.cos(phi0) * np.cos(phi) * np.cos(lam)
    k = np.sqrt(2.0 / np.maximum(1.0 + cosc, 1e-12))
    x = R_KM * k * np.cos(phi) * np.sin(lam)
    y = R_KM * k * (math.cos(phi0) * np.sin(phi) - math.sin(phi0) * np.cos(phi) * np.cos(lam))
    return x, y, cosc


def laea_inv(x, y, lon0, lat0):
    rho = math.hypot(x, y)
    if rho < 1e-9:
        return lon0, lat0
    c = 2 * math.asin(min(1.0, rho / (2 * R_KM)))
    phi0 = math.radians(lat0)
    lat = math.asin(math.cos(c) * math.sin(phi0) + y * math.sin(c) * math.cos(phi0) / rho)
    lon = math.radians(lon0) + math.atan2(
        x * math.sin(c), rho * math.cos(phi0) * math.cos(c) - y * math.sin(phi0) * math.sin(c))
    return (math.degrees(lon) + 540) % 360 - 180, math.degrees(lat)


def project_ring(coords, lon0, lat0):
    a = np.asarray(coords)[:, :2]
    x, y, cosc = laea(a[:, 0], a[:, 1], lon0, lat0)
    return np.column_stack([x, y]), cosc.min()


def project_polygon(poly, lon0, lat0):
    """Project one lon/lat polygon; None if it wraps round the antipode."""
    ext, mn = project_ring(poly.exterior.coords, lon0, lat0)
    if mn < -0.985:                              # a vertex within ~10° of antipode
        return None
    holes = [project_ring(r.coords, lon0, lat0)[0] for r in poly.interiors]
    return Polygon(ext, holes)


def polys(geom):
    if geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return list(geom.geoms)
    return [g for g in getattr(geom, "geoms", []) if isinstance(g, Polygon)]


def lines(geom):
    if geom.is_empty:
        return []
    if isinstance(geom, LineString):
        return [geom]
    if isinstance(geom, MultiLineString):
        return list(geom.geoms)
    return [g for gg in getattr(geom, "geoms", []) for g in lines(gg)]


def antipode_inside(poly, lon0, lat0):
    alon = (lon0 + 360) % 360 - 180
    return poly.contains(Point(alon, -lat0))


# ----------------------------------------------------------------- svg paths
def fmt(v):
    s = f"{v:.1f}"
    return s[:-2] if s.endswith(".0") else s


def ring_path(coords):
    pts = list(coords)
    if len(pts) > 1 and pts[0] == pts[-1]:
        pts = pts[:-1]
    if len(pts) < 3:
        return ""
    out, last = [], None
    for x, y in pts:
        p = (fmt(x), fmt(y))
        if p != last:
            out.append(p)
            last = p
    if len(out) < 3:
        return ""
    return "M" + "L".join(f"{a} {b}" for a, b in out) + "Z"


def line_path(coords):
    out, last = [], None
    for x, y in coords:
        p = (fmt(x), fmt(y))
        if p != last:
            out.append(p)
            last = p
    if len(out) < 2:
        return ""
    return "M" + "L".join(f"{a} {b}" for a, b in out)


def poly_path(geoms):
    parts = []
    for g in geoms:
        for p in polys(g):
            parts.append(ring_path(p.exterior.coords))
            parts.extend(ring_path(r.coords) for r in p.interiors)
    return "".join(parts)


# ----------------------------------------------------------------- load data
def load(data_dir):
    fc = json.loads((data_dir / "ne_50m_admin_0_countries.geojson").read_text())
    feats = {}
    for f in fc["features"]:
        p = f["properties"]
        code = p["ISO_A3_EH"] if p["ISO_A3_EH"] not in ("-99", None) else p["ADM0_A3"]
        geom = shapely.make_valid(shape(f["geometry"]))
        feats[p["ADM0_A3"]] = dict(code=code, a3=p["ADM0_A3"], props=p, geom=geom)
    for src, dst in MERGE_INTO.items():
        feats[dst]["geom"] = unary_union([feats[dst]["geom"], feats.pop(src)["geom"]])
    lc = json.loads((data_dir / "ne_50m_admin_0_boundary_lines_land.geojson").read_text())
    blines = []
    for f in lc["features"]:
        dashed = f["properties"].get("FEATURECLA") in DASHED_CLASSES
        blines.extend((ln, dashed) for ln in lines(shape(f["geometry"])))
    return feats, blines


# --------------------------------------------------------------- view frame
def frame_for(geom):
    """Pick projection centre and km bounding box for a country's main body."""
    parts = polys(geom)

    def approx_area(p):
        return p.area * math.cos(math.radians(p.centroid.y))

    main = max(parts, key=approx_area)
    lon0, lat0 = main.centroid.x, main.centroid.y
    for _ in range(2):
        proj = [project_polygon(p, lon0, lat0) for p in parts]
        proj = [p for p in proj if p is not None]
        pmain = max(proj, key=lambda p: p.area)
        minx, miny, maxx, maxy = pmain.bounds
        diag = math.hypot(maxx - minx, maxy - miny)
        thresh = max(700.0, min(0.5 * diag, 1000.0))
        # single-linkage cluster of nearby parts (keeps island chains together)
        cluster, rest = [pmain], [p for p in proj if p is not pmain]
        grew = True
        while grew:
            grew = False
            for p in list(rest):
                d = min(p.distance(c) for c in cluster)
                if d < thresh and not (p.area < 0.01 * pmain.area and d > 250):
                    cluster.append(p)
                    rest.remove(p)
                    grew = True
        minx, miny, maxx, maxy = MultiPolygon(cluster).bounds
        lon0, lat0 = laea_inv((minx + maxx) / 2, (miny + maxy) / 2, lon0, lat0)
    proj = [p for p in (project_polygon(p, lon0, lat0) for p in parts) if p is not None]
    total_area = sum(p.area for p in proj)
    return lon0, lat0, (minx, miny, maxx, maxy), total_area


def build_view(target_key, feats, blines):
    tgt = feats[target_key]
    lon0, lat0, (minx, miny, maxx, maxy), area = frame_for(tgt["geom"])
    w, h = maxx - minx, maxy - miny
    pad = 0.18
    # small countries get a wider frame so their neighbours are visible
    w, h = max(w * (1 + 2 * pad), 950.0), max(h * (1 + 2 * pad), 710.0)
    cx, cy = (minx + maxx) / 2, (miny + maxy) / 2
    aspect = VIEW_W / VIEW_H
    if w / h < aspect:
        w = h * aspect
    else:
        h = w / aspect
    scale = VIEW_W / w                     # px per km

    def to_px(g):
        return shapely.transform(g, lambda a: np.column_stack(
            [(a[:, 0] - cx) * scale + VIEW_W / 2, (cy - a[:, 1]) * scale + VIEW_H / 2]))

    margin = 12
    clip = box(-margin, -margin, VIEW_W + margin, VIEW_H + margin)
    kmclip = box(cx - (w / 2 + margin / scale), cy - (h / 2 + margin / scale),
                 cx + (w / 2 + margin / scale), cy + (h / 2 + margin / scale))

    def render_geom(geom):
        out = []
        for p in polys(geom):
            if antipode_inside(p, lon0, lat0):
                continue
            pp = project_polygon(p, lon0, lat0)
            if pp is None or not pp.intersects(kmclip):
                continue
            if not pp.is_valid:
                pp = shapely.make_valid(pp)
            g = to_px(pp.intersection(kmclip)).intersection(clip)
            g = g.simplify(SIMPLIFY_PX, preserve_topology=True)
            if not g.is_empty:
                out.append(g)
        return out

    land = []
    for key, f in feats.items():
        if key != target_key:
            land.extend(render_geom(f["geom"]))
    target = render_geom(tgt["geom"])

    solid, dashed = [], []
    for ln, is_dashed in blines:
        a = np.asarray(ln.coords)
        x, y, cosc = laea(a[:, 0], a[:, 1], lon0, lat0)
        if cosc.min() < 0:
            continue
        pl = LineString(np.column_stack([x, y]))
        if not pl.intersects(kmclip):
            continue
        g = to_px(pl.intersection(kmclip)).simplify(SIMPLIFY_PX)
        for seg in lines(g):
            (dashed if is_dashed else solid).append(line_path(seg.coords))

    return dict(
        land=poly_path(land),
        borders="".join(solid),
        disputed="".join(dashed),
        target=poly_path(target),
        km_per_px=round(1 / scale, 3),
    ), (lon0, lat0), area


# ----------------------------------------------------------------- distractors
def gc_km(a, b):
    lon1, lat1 = map(math.radians, a)
    lon2, lat2 = map(math.radians, b)
    h = math.sin((lat2 - lat1) / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin((lon2 - lon1) / 2) ** 2
    return 2 * R_KM * math.asin(min(1.0, math.sqrt(h)))


def rank_distractors(ids, info, neighbours, keep=8):
    out = {}
    for i in ids:
        scored = []
        for j in ids:
            if j == i:
                continue
            d = gc_km(info[i]["center"], info[j]["center"])
            size = abs(math.log(info[j]["area"] / info[i]["area"]))
            s = d * (1 + 0.25 * size)
            if j in neighbours[i]:
                s *= 0.25
            scored.append((s, j))
        scored.sort()
        out[i] = [j for _, j in scored[:keep]]
    return out


def main():
    data_dir = Path(sys.argv[1])
    out_html = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parent.parent / "country-quiz.html"
    feats, blines = load(data_dir)

    eligible = {}
    for key, f in feats.items():
        code = f["code"]
        if f["props"]["TYPE"] != "Dependency" and (code in UN_MEMBERS or code in UN_OBSERVERS):
            eligible[key] = f
    found = {f["code"] for f in eligible.values()}
    missing = (UN_MEMBERS | UN_OBSERVERS) - found
    if missing:
        print("WARNING: not in Natural Earth 50m:", sorted(missing), file=sys.stderr)

    countries, info, skipped = [], {}, []
    for key in sorted(eligible, key=lambda k: eligible[k]["props"]["NAME"]):
        f = eligible[key]
        view, center, area = build_view(key, feats, blines)
        name = NAME_OVERRIDES.get(f["code"], f["props"]["NAME_EN"] or f["props"]["NAME"])
        if area < MIN_AREA_KM2:
            skipped.append((name, round(area)))
            continue
        region = REGION_OF_CONTINENT.get(f["props"]["CONTINENT"])
        idx = len(countries)
        countries.append(dict(name=name, code=f["code"], region=region, view=view))
        info[idx] = dict(key=key, center=center, area=area, region=region, code=f["code"])
        print(f"{name:40s} {region or '-':9s} {area:12,.0f} km²  {len(json.dumps(view)) / 1024:6.1f} KB",
              file=sys.stderr)

    # land neighbours (shared border) among quiz countries
    geoms = {i: feats[info[i]["key"]]["geom"] for i in info}
    tree_ids = list(geoms)
    tree = shapely.STRtree([geoms[i] for i in tree_ids])
    neighbours = {i: set() for i in info}
    for i in info:
        for k in tree.query(geoms[i].buffer(0.02)):
            j = tree_ids[k]
            if j != i and geoms[i].distance(geoms[j]) < 0.02:
                neighbours[i].add(j)

    modes = []
    for mode, label in MODES:
        ids = [i for i in info if (mode == "world" or info[i]["region"] == mode)
               and (mode == "world" or info[i]["code"] not in UN_OBSERVERS)]
        modes.append(dict(id=mode, label=label, ids=ids,
                          distractors=rank_distractors(ids, info, neighbours)))

    data = dict(
        source="Natural Earth 1:50m Admin 0 – Countries (public domain), naturalearthdata.com",
        view=[VIEW_W, VIEW_H],
        countries=[dict(n=c["name"], v=c["view"]) for c in countries],
        modes=modes,
    )
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    template = (Path(__file__).resolve().parent / "quiz_template.html").read_text()
    out_html.write_text(template.replace("/*__QUIZ_DATA__*/null", payload))
    print(f"\nLeft out as too small (<{MIN_AREA_KM2:,} km²): " +
          ", ".join(f"{n} ({a:,})" for n, a in skipped), file=sys.stderr)
    for m in modes:
        print(f"{m['label']}: {len(m['ids'])} countries", file=sys.stderr)
    print(f"Wrote {out_html} ({out_html.stat().st_size / 1024:,.0f} KB)", file=sys.stderr)


if __name__ == "__main__":
    main()
