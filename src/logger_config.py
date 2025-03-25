import logging

# Configure logging (set up once)
logging.basicConfig(filename="helical_fmo.log",
                    filemode="a",
                    format="%(asctime)s - %(levelname)s - %(message)s",
                    level=logging.DEBUG)



def get_logger(name):
    return logging.getLogger(name)
