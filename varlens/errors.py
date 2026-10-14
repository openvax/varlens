"""Actionable input errors shared by the Python and command-line APIs."""

from functools import wraps


class InputError(ValueError):
    """The input cannot be processed without losing or misinterpreting data."""


def input_errors(parser):
    """Give CLI users a concise diagnostic, while Python APIs raise normally."""
    def decorate(function):
        @wraps(function)
        def run(*args, **kwargs):
            try:
                return function(*args, **kwargs)
            except InputError as error:
                parser.error(str(error))
        return run
    return decorate
