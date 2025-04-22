from abc import ABC, abstractmethod
from typing import Dict, Union

class Controller(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def run_controller(self):
        pass

    @abstractmethod
    def validate_controller(self, config_section: Dict[str, Union[str, int, float, bool]]):
        pass

