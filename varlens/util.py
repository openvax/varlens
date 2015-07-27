try:
    # Python 2
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

try:
    # Python 2
    from urlparse import parse_qsl
except ImportError:
    from urllib.parse import parse_qsl

def parse_url_fragment(url):
    parsed = urlparse(url)
    try:
        if parsed.fragment:
            # If our fragment begins with an '&' symbol, we ignore it. 
            unparsed_fragment = parsed.fragment
            if unparsed_fragment.startswith("&"):
                unparsed_fragment = unparsed_fragment[1:]
            fragment = parse_qsl(unparsed_fragment, strict_parsing=True)
        else:
            fragment = []
    except ValueError as e:
        raise ValueError("Couldn't parse fragment '%s': %s" % (
            parsed.fragment, e))

    return (parsed._replace(fragment='').geturl(), fragment)

def string_to_boolean(s):
    value = str(s).lower()
    if value in ("true", "1"):
        return True
    elif value in ("false", 0):
        return False
    raise ValueError("Not a boolean string: %s" % s)