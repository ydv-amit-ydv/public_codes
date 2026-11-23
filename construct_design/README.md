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

## Two directories

```
construct_design/
  app.py
  gene-construct-builder-server.html
  gene-construct-builder.html
  data/          <- constructs people save (git-ignored)
  library/       <- shared constructs everyone can read (committed)
```

Each construct is one plain, pretty-printed `.json` file named after it, e.g. `my-secreted-fusion.json` — copy, back up, version-control, hand-edit or delete them freely. Drop a file into either folder and it appears in the list (the page re-reads both every 15 seconds). Renaming a construct in the UI renames its file.

Change either with `--data <dir>` / `--library <dir>`.

## Shared library

`library/` is a common, read-only set everyone on the LAN can see, in its own section of the page:

- **Open** loads it into the editor as an unsaved construct, so your first Save writes a copy into `data/` and the library file is untouched.
- **Copy to my constructs** copies it into `data/` without opening it.
- **Download** fetches the raw file.

The server never writes to `library/`; `PUT` and `DELETE` on a library path return 405. To add a shared construct, put a `.json` file in that folder (or save one in `data/` and move it over). Four starters ship with it: a generic mRNA, a tPA-leader secreted fusion, a bicistronic IRES cassette and an exon/intron gene.

## API

| Method | Path | Does |
|---|---|---|
| GET | `/api/constructs` | list all, with the data directory path |
| GET | `/api/constructs/<id>` | read one (`id` = filename without `.json`) |
| POST | `/api/constructs` | create from `{name, settings, blocks}` |
| PUT | `/api/constructs/<id>` | update, renaming the file if the name changed |
| DELETE | `/api/constructs/<id>` | delete |
| GET | `/api/library` | list the shared read-only library |
| GET | `/api/library/<id>` | read one shared construct |
| GET | `/files/<id>.json` | download the raw file |
| GET | `/files/library/<id>.json` | download a shared one |

Every read endpoint is open to anyone who can reach the port, so other tools and scripts on the LAN can pull a construct straight out of either folder.

## Security

There is **no login**: anyone who can reach the port can read, change and delete every construct in the folder. Run it only on a trusted LAN and don't port-forward it. Requests cannot escape the two directories, and the server only reads and writes `.json` files inside them — `library/` is read-only even through the API. Keep a backup — Delete is immediate and permanent.

Flask's built-in server is fine for a few people on a LAN. For heavier use: `pip install waitress && waitress-serve --host 0.0.0.0 --port 8080 app:app` (with waitress the directories stay at their defaults, `./data` and `./library`).
