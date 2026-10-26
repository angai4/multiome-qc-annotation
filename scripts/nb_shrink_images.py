"""Git clean filter: shrink the plots embedded in a notebook so the committed copy stays small enough for GitHub.

The notebooks save plots at 300 dpi (up to 14,000 px wide, ~190 MB per notebook). Each image/png output wider or
taller than MAX_SIDE is resized to fit MAX_SIDE and stored as PNG or JPEG (quality 90), whichever is smaller; dense
scatter plots compress far better as JPEG, line plots as PNG. Everything else is copied unchanged. Local notebooks
keep their full-resolution images; only what git stores is shrunk. Needs Pillow (the multiome env).

  git config filter.nbshrink.clean "<multiome env python> scripts/nb_shrink_images.py"
  git config filter.nbshrink.smudge cat
"""

import base64
import io
import json
import sys

from PIL import Image

MAX_SIDE = 1600
JPEG_QUALITY = 90
Image.MAX_IMAGE_PIXELS = None  # the notebooks' own plots, up to ~110 Mpx


def shrink(png_b64):
    """(mime type, base64 data) of the smaller of PNG and JPEG, or None if the image is already small."""
    with Image.open(io.BytesIO(base64.b64decode(png_b64))) as im:
        if max(im.size) <= MAX_SIDE:
            return None
        im = im.copy()
    im.thumbnail((MAX_SIDE, MAX_SIDE), Image.LANCZOS)
    if im.mode != "RGB":  # figures are drawn on white, so flattening the alpha channel onto white changes nothing
        rgb = Image.new("RGB", im.size, "white")
        rgb.paste(im, mask=im.convert("RGBA").getchannel("A"))
        im = rgb
    png, jpeg = io.BytesIO(), io.BytesIO()
    im.save(png, "PNG", optimize=True)
    im.save(jpeg, "JPEG", quality=JPEG_QUALITY, optimize=True)
    mime, data = min(("image/png", png.getvalue()), ("image/jpeg", jpeg.getvalue()), key=lambda x: len(x[1]))
    return mime, base64.b64encode(data).decode()


nb = json.load(sys.stdin)
for cell in nb["cells"]:
    for output in cell.get("outputs", []):
        data = output.get("data", {})
        if "image/png" not in data:
            continue
        shrunk = shrink(data["image/png"])
        if shrunk is None:
            continue
        mime, b64 = shrunk
        del data["image/png"]
        data[mime] = b64
        meta = output.get("metadata", {})
        if "image/png" in meta:
            meta[mime] = meta.pop("image/png")

# nbformat's own layout, so diffs against Jupyter-saved files stay small
sys.stdout.write(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
