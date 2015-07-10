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