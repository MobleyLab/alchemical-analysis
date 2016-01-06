from logger import logger



class GeneralException(Exception):
    pass

class ParserException(Exception):
    pass


def log_and_raise(msg, excpt=ParserException):
    """Log and raise and exception."""

    logger.error(msg)
    raise excpt('ERROR: ' + msg)
