"""_summary_

This script contains the definition of Cell, CancerCell, Hepatocyte classes.
Main methods of the classes are for getting and setting attributes.
    
"""

from typing import Dict

class Cell:
    def __init__(self, cell_attributes: Dict):
        self.attributes = cell_attributes 
        
    def get_attributes(self) -> Dict:
        return self.attributes
    
    def set_attributes(self, updated_attributes: Dict):
        self.attributes = updated_attributes
        
class CancerCell(Cell):
    def __init__(self, cell_attributes: Dict):
        super().__init__(cell_attributes)

class Hepatocyte(Cell):
    def __init__(self, cell_attributes: Dict):
        super().__init__(cell_attributes)