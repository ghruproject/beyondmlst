"""Stable location colours shared by network and tree presentations."""


def country_palette(countries):
    """Stable country colours independent of the selected view."""
    import colorsys
    import hashlib

    return {
        country: "#"
        + "".join(
            f"{round(v * 255):02x}"
            for v in colorsys.hsv_to_rgb(
                int.from_bytes(hashlib.sha256(country.encode()).digest()[:4], "big") / 2**32,
                0.58,
                0.72,
            )
        )
        for country in countries
    }
