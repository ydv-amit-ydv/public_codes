# Home Designer

Our own floor-plan and interior design tool. You draw walls and rooms in 2D, furnish them from a
**337-item catalog**, annotate places with notes, and see the result in a realistic **3D view**.
Designs are saved as JSON files in [`designs/`](designs/), so they live in git with everything else.

It runs locally in the browser. The only requirement is Python 3; there's no build step and no npm install.

## Start

```bash
cd home_designer && bash start.sh
```

This opens <http://127.0.0.1:8765>. On Windows, double-click `start.bat`. Other options:

```bash
python3 server.py --port 9000                # different port
python3 server.py --host 0.0.0.0             # open from a phone/tablet on the same Wi-Fi
HOME_DESIGNER_DATA=~/my-designs python3 server.py   # keep designs somewhere else
```

When you first run it, the sample design (`designs/sample_2bhk.json`) opens: a 2BHK house with a garden,
parking and a terrace.

## What it does

**Drawing**
- **Walls** (W): click to place each corner. Shift snaps to 45°, Alt turns snapping off.
- **Exact lengths, fast:** while drawing, type the length (`4.5`, `450cm`, `12'6`, `15ft`) and press Enter.
  To resize an existing wall, select it, type the new length and press Enter. The joined walls move with it, so the room just gets bigger.
  With the room tool, type `4x3.5` and press Enter.
- **Two wall types:** *Exterior* (default 23 cm / 9″) and *Interior* (default 11.5 cm / 4.5″). Switch with the
  Ext/Int button or press `E`. Change either thickness in Design settings and every wall of that type updates. *Auto-detect* marks
  walls that touch the outside as exterior.
- **Rooms:** room + walls (B), area only (A) for gardens, parking and so on, or a polygon (P) for L-shaped rooms.
- **Doors and windows** (18 types) snap into walls and slide along them. You can flip the swing and the hinge.
- **Multiple floors**, with the floor below drawn faintly as a guide. Undo/redo, duplicate, nudge, and metres or feet.

**Annotating places**
- Each room has a name, type (44 presets), floor finish, wall paint and notes. Area and perimeter are calculated automatically.
- 📍 **Note pins** (numbered), text labels and dimension lines.
- Every item, door, window and wall has a *label* and *notes* field (brand, price, link, status…).
- **Places & notes** tab: a searchable list of rooms, notes and annotated items.
- **Schedule** tab: item counts by category. Exports to CSV, or a printable/PDF sheet per floor.

**Catalog: 337 items in 21 categories.** Doors & windows, structure (columns, stairs, ramps, railings), living,
dining, kitchen (counters, hob, chimney, RO, gas cylinder…), bedroom, bathroom, laundry & utility, study,
kids, pooja (mandir, diya stand), storage, lighting & fans, electrical & HVAC (switchboards, sockets, split AC,
DB, CCTV, inverter, EV charger), plumbing (tanks, sump, RWH pit), decor, outdoor & garden (trees, tulsi,
swing, pergola, pool, gate, compound wall), vehicles, fitness & entertainment, pets.
Items marked as wall items (TV, AC, wardrobes, cabinets…) snap flat against the nearest wall.

**3D view**
- three.js with PBR materials, image-based lighting, GTAO ambient occlusion, ACES tone mapping and soft
  shadows. The sun direction is adjustable.
- **Real photo textures** (CC0, [Poly Haven](https://polyhaven.com)) in [`web/textures/`](web/textures/): wood, parquet,
  marble, tiles, terrazzo, terracotta, granite/kota, carpet, grass, pavers, concrete, plus fabric, leather, wood
  veneer and plaster for furniture and walls. Textures are mapped at their real-world size.
  Each room's floor finish picks the texture automatically, or you can choose one under *3D floor*.
- **Real 3D models:** select an item, then *3D model → Browse free models* to choose from **522 CC0 models** from Poly Haven
  (sofas, chairs, lamps, plants, tables…), or *Upload .glb* for your own (from Sketchfab, Kenney, Blender…).
  Models download into `designs/models/` (commit that folder with your designs). *Use for every …* applies a
  model to all items of that type.
- **Walk** mode (first person, WASD), *low walls* cut-away, top view, and HD screenshots (📷).

**Trace a sketch, photo or existing plan**
- *Import / Export → Trace over an image*, or *Design settings → Trace image*. Upload a file or paste a web URL.
- Click both ends of something whose length you know and type that length. The image is scaled to match.
- Lower the opacity and draw walls over it. The image is stored inside the design.
- You can also send Claude a photo or sketch and ask it to write the design JSON directly. The format is described below.

## Saving

| How | Where |
|---|---|
| **Save** (Ctrl+S) | `designs/<name>.json`. Commit these files to keep them. |
| Every save | The previous version is kept in `designs/.history/` (last 20 versions, git-ignored). |
| Auto-draft | Your current work is also kept in the browser and restored when you reopen the page. |
| No server (static hosting) | Designs go to browser storage. Use *Export design* to download them. |

## Keys

`V` select · `W` wall · `B` room · `A` area · `P` polygon · `N` note · `T` label · `D` dimension · `H` pan ·
`E` exterior/interior walls · type a number + `Enter` for an exact length · `R` rotate (Shift+R: 15°) ·
`Del` delete · `Ctrl+D` duplicate · `Ctrl+Z`/`Ctrl+Y` undo/redo · arrows nudge 1 cm (Shift: 10 cm) ·
`F` fit · `Space`+drag pan · `1`/`2`/`3` 2D / split / 3D · `Esc` cancel.

## Design file format

Plain JSON. All lengths are in metres. The plan's x axis points right and y points down. Angles are in degrees.

```jsonc
{
  "format": "home-designer", "version": 1, "name": "My home",
  "settings": { "units": "m", "grid": 0.1, "extThickness": 0.23, "intThickness": 0.115, "wallHeight": 3 },
  "floors": [{
    "name": "Ground floor", "elevation": 0, "height": 3,
    "walls":    [{ "id": "w1", "x1": 0, "y1": 0, "x2": 6, "y2": 0, "thickness": 0.23, "kind": "exterior", "height": 3 }],
    "openings": [{ "id": "o1", "wallId": "w1", "catalogId": "window_std", "t": 0.5, "width": 1.2, "height": 1.2, "sill": 0.9, "swing": 1 }],
    "rooms":    [{ "id": "r1", "name": "Living", "type": "Living Room", "points": [[0,0],[6,0],[6,4],[0,4]], "color": "#fff3e0", "floorFinish": "vitrified tiles", "notes": "" }],
    "items":    [{ "id": "i1", "catalogId": "sofa3", "x": 3, "y": 3.4, "rot": 180, "w": 2.1, "d": 0.9, "h": 0.85, "elev": 0, "label": "", "notes": "", "model3d": "models/polyhaven/sofa_02/sofa_02_1k.gltf" }],
    "notes":    [{ "id": "n1", "kind": "pin", "x": 1, "y": 1, "title": "Socket", "text": "Add 16A socket" }],
    "underlay": { "src": "asset:img_x", "x": 0, "y": 0, "mpp": 0.01, "opacity": 0.5 }
  }],
  "assets": { "img_x": "data:image/jpeg;base64,…" }
}
```

An opening's `t` is its position along the wall (0 = start, 1 = end). An item's `x`/`y` is its centre, and its front faces +y
when `rot` is 0. Catalog ids are listed in [`web/js/catalog.js`](web/js/catalog.js).

## Files

```
home_designer/
  start.sh / start.bat     launch
  server.py                stdlib server: static files, designs API, image fetch, Poly Haven proxy
  designs/                 saved designs (*.json) + models/ (downloaded/uploaded 3D models)
  web/index.html           the app
  web/js/catalog.js        337 items + room presets
  web/js/plan.js           2D editor
  web/js/view3d.js         3D renderer
  web/js/textures.js       photo + procedural textures
  web/js/main.js           panels, tools, save/load
  web/textures/            CC0 photo textures (Poly Haven) + manifest.json
  web/js/vendor/           three.js r160 (MIT) + addons
```

## Credits and licences

- [three.js](https://threejs.org) r160, MIT licence (`web/js/vendor/THREE_LICENSE.txt`).
- Photo textures and downloadable models from [Poly Haven](https://polyhaven.com), CC0. They're free for any use and need no attribution,
  but authors are listed in `web/textures/manifest.json`.
