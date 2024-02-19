import logging

from opentps.core.data.plan import RTPlan

logger = logging.getLogger(__name__)

class PhotonPlan(RTPlan):
    """
            Class for storing the data of a single PhotonPlan. Inherits from RTPlan.

            Attributes
            ----------

    """
    def __init__(self, name="PhotonPlan", patient=None):
        super().__init__(name=name, patient=patient)