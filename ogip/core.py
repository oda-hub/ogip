import logging

logger=logging.getLogger()

import ogip.spec

def open_something(fn, allow_many=False):
    options = []

    for c in ogip.spec.PHAI, ogip.spec.RMF:
        try:
            options.append(c.from_file_name(fn))
            logger.info(f"opened {fn} as {options[-1]}")
        except Exception as e:
            logger.debug(e)

    return options

