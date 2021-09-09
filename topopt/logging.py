import logging
logger = logging.getLogger('topopt')

# Logger Class copied from https://github.com/ouyang-w-19/decogo
class TopOptLogger:
    """Class which is responsible for the logging. It sets up
    the logs
    """

    def __init__(self, level, file_name=None):
        """Constructor method"""
        if level == 'debug':
            level = logging.DEBUG
        elif level == 'info':
            level = logging.INFO

        if file_name is None:
            self._set_up_into_console(level)
        else:
            self._set_up_into_file(file_name, level)

    def _set_up_into_file(self, file_name, level):
        """Sets up the logger into file

        :param file_name: Name of the file
        :type file_name: str
        :param level: logger level
        :type level: logging.DEBUG, logging.INFO etc
        """
        logger.setLevel(level)

        ch = logging.FileHandler(file_name, mode='w')
        ch.setLevel(logging.DEBUG)

        formatter = logging.Formatter('%(message)s')

        ch.setFormatter(formatter)
        logger.addHandler(ch)

    def _set_up_into_console(self, level):
        """Sets up the logger into console

        :param level: logger level
        :type level: logging.DEBUG, logging.INFO etc
        """

        logger.setLevel(level)

        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        formatter = logging.Formatter('%(asctime)s [%(levelname)-5s] %(message)s')

        ch.setFormatter(formatter)
        logger.addHandler(ch)

    def clean_up(self):
        """Cleans all logger handlers"""
        for i, handler in enumerate(logger.handlers):
            logger.handlers[i].close()
            logger.removeHandler(handler)
