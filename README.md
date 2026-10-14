# Name That Country

A map quiz for 7th grade world geography. Open **[`country-quiz.html`](country-quiz.html)** in any browser (Chrome on a Chromebook works). It's one self-contained file with no login, no internet needed after it loads, and no student data collected or saved.

- Pick a region: Africa, Asia, Europe, North & Central America, South America, or the Whole World.
- Each question shows one country highlighted in orange on a map, with neighboring countries in gray for context. Students choose from four names.
- Wrong answers are plausible: mostly bordering countries, then nearby countries of similar size from the same region.
- Students get instant feedback (✓/✗). If they miss, the page shows the correct answer before moving on.
- A round is 10 questions with no repeats. At the end students see their score, a review card for each country they missed, and a **Play Again** button.
- Keyboard: keys **1–4** answer, **Enter** goes to the next question.

## Map data and choices

- Boundaries come from [Natural Earth](https://www.naturalearthdata.com/) 1:50m Admin 0 Countries (public domain). Each country is drawn in an equal-area projection centered on it, so its shape isn't stretched.
- Countries quizzed: the 193 UN member states, plus Palestine (a UN observer state) in Whole World mode only, minus those smaller than 5,000 km² (e.g., Seychelles, Comoros, Cabo Verde, São Tomé and Príncipe, Mauritius, Singapore, Luxembourg, Caribbean and Pacific microstates). That leaves 165 countries; Oceania (6 countries left) appears in Whole World mode.
- Somaliland is shown as part of Somalia, and Northern Cyprus as part of Cyprus.
- Western Sahara, Kosovo and Taiwan are drawn as gray areas: never highlighted and never quizzed. Disputed borders and lines of control (e.g., Morocco–Western Sahara, Kashmir) are dashed.

## Rebuilding

`country-quiz.html` is generated. To change the template, names, region groupings or size cutoff:

```sh
pip install shapely numpy
B=https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson
mkdir -p ne && for f in ne_50m_admin_0_countries ne_50m_admin_0_boundary_lines_land; do curl -sSo ne/$f.geojson $B/$f.geojson; done
python3 tools/build_quiz.py ne
```
