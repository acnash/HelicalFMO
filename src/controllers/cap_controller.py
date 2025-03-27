from src.controllers.controller import Controller


class CapController(Controller):

    def __init__(self, input_file, input_directory):
        self.input_file = input_file
        self.input_directory = input_directory

        self.universe_list = []

    def validate_inputs(self):
        if self.input_file:
            pass
        else:
            pass

    def run_controller(self):
        pass

