# Pipeline figure

Source and exported versions of the pipeline schematic used in the README and as
a supplementary methods figure.

| File | What it is |
|---|---|
| `pipeline_figure.html` | Self-contained source (HTML/CSS). Edit this, then re-export. |
| `pipeline.png` | Raster export used in the top-level `README.md`. |
| `pipeline.pdf` | Vector export for editing in Adobe Illustrator (live text, embedded fonts). |

## Regenerating the exports

The exports were produced with Google Chrome in headless mode (macOS path shown;
adjust for your platform) and, for the SVG, the Mermaid CLI.

PNG (as used in the README):

```bash
CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
"$CHROME" --headless=new --disable-gpu --hide-scrollbars \
  --force-device-scale-factor=2 --window-size=1300,5145 \
  --screenshot=pipeline.png "file://$PWD/pipeline_figure.html"
```

PDF (vector, for Illustrator):

```bash
"$CHROME" --headless=new --disable-gpu --no-pdf-header-footer \
  --print-to-pdf=pipeline.pdf "file://$PWD/pipeline_figure.html"
```

The window height must match the rendered document, or the export gains a band
of blank paper at the bottom (or clips the last line). After editing the figure,
re-measure it: screenshot at an over-generous height, find the last row that is
not page background, halve it (the export is at 2x) and add the 24px bottom
padding. The current figure is 5145 CSS px tall.

The single-page size and print colours are set by the `@media print` block near
the end of the `<style>` in `pipeline_figure.html`.

Opening `pipeline.pdf` in Illustrator: use File > Open. Text is editable; the
export nests objects in clip groups, so ungroup (Cmd+Shift+G) once or twice to
reach individual boxes.
