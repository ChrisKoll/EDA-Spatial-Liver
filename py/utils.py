# Standard imports
import logging
import os
import sys


def verify_path(path: str, is_dir: bool = False) -> str:
    """
    Verify that the given path exists and is of the expected type (file or directory).

    :param path: The path to verify.
    :type path: str
    :param is_dir: If True, checks that the path is a directory. If False, checks that the path is a file.
    :type is_dir: bool, optional

    :returns: The verified path if it exists and matches the expected type.
    :rtype: str

    :raises FileNotFoundError: If the path does not exist or does not match the expected type.
    """

    if is_dir:
        if not os.path.isdir(path):  # Verify if path truly leads to directory
            err_msg = f"""
            Path must be a directory!\n
            `{path}` leads to a file.
            """
            raise FileNotFoundError(err_msg)
        else:
            return path
    else:
        if not os.path.isfile(path):  # Verify if path truly leads to file
            err_msg = f"""
            Path must lead to file!\n 
            `{path}` leads to a directory.
            """

            raise FileNotFoundError(err_msg)
        else:
            return path


def setup_logging(path_to_log: str) -> logging.Logger:
    """
    Set up and configure logging to both console and file.

    :param path_to_log: Path to the log file where logs will be written.
    :type path_to_log: str

    :returns: Configured logger instance.
    :rtype: logging.Logger
    """

    logger = logging.getLogger("data_integration")
    logger.setLevel(logging.DEBUG)  # Lowest log level (logs everything)

    # Custom formatter
    formatter = logging.Formatter(
        "%(asctime)s || LEVEL: %(levelname)s |> %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)  # Set handler specific log level
    console_handler.setFormatter(formatter)  # Add custom formatter to handler
    logger.addHandler(console_handler)  # Add handler to the logger

    # Add file handler
    file_handler = logging.FileHandler(path_to_log)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger
