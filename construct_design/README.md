# Gene Construct Builder

| File | Use |
|---|---|
| `gene-construct-builder.html` | Open straight in a browser. Constructs are saved in that browser only. |
| `app.py` + `gene-construct-builder-server.html` | Flask server: run on one machine, everyone on the LAN shares the same constructs. |

## Run the server

```bash
pip install flask
```

```bash
python3 app.py
```

It prints the address to share, e.g. `http://192.168.1.40:8080/`.

Options: `--data <dir>` (default: the folder `app.py` is in), `--port <n>` (default 8080), `--host <addr>` (default `0.0.0.0`, the whole LAN).

## One directory

`app.py`, `gene-construct-builder-server.html` and the constructs all live in the same folder. Each construct is one plain, pretty-printed `.json` file named after it, e.g. `bicistronic-tpa-test.json` — copy, back up, version-control, hand-edit or delete them freely. Drop a file in from elsewhere and it appears in the list (the page re-reads the folder every 15 seconds). Renaming a construct in the UI renames its file.

To keep constructs apart from the code: `python3 app.py --data ./constructs`

## API

| Method | Path | Does |
|---|---|---|
| GET | `/api/constructs` | list all, with the data directory path |
| GET | `/api/constructs/<id>` | read one (`id` = filename without `.json`) |
| POST | `/api/constructs` | create from `{name, settings, blocks}` |
| PUT | `/api/constructs/<id>` | update, renaming the file if the name changed |
| DELETE | `/api/constructs/<id>` | delete |
| GET | `/files/<id>.json` | download the raw file |

## Security

There is **no login**: anyone who can reach the port can read, change and delete every construct in the folder. Run it only on a trusted LAN and don't port-forward it. Requests cannot escape the data directory, and the server only reads and writes `.json` files inside it. Keep a backup — Delete is immediate and permanent.

Flask's built-in server is fine for a few people on a LAN. For heavier use: `pip install waitress && waitress-serve --host 0.0.0.0 --port 8080 app:app` (set the data directory with `--data` only via `app.py`; with waitress it defaults to the app folder).
