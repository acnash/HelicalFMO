from abc import ABC, abstractmethod
class Controller(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def run_controller(self):
        pass

