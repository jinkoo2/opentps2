from PyQt5.QtWidgets import QWidget, QVBoxLayout, QComboBox, QLabel, QPushButton

from opentps.core.data.images import DoseImage, CTImage
from opentps.core.data.images._image3D import Image3D
from opentps.core.data._patient import Patient


def _dose_like_data(patient):
    """DoseImage and Image3D that are not CT (so MHD-loaded doses appear)."""
    if patient is None:
        return []
    doses = patient.getPatientDataOfType(DoseImage)
    others = [
        im for im in patient.images
        if isinstance(im, Image3D) and not isinstance(im, (CTImage, DoseImage))
    ]
    return list(doses) + others


def _to_dose_image(data):
    """Use as DoseImage; convert Image3D if needed."""
    if data is None:
        return None
    if isinstance(data, DoseImage):
        return data
    return DoseImage.fromImage3D(data)


class DoseComparisonPanel(QWidget):
    def __init__(self, viewController):
        QWidget.__init__(self)

        self._patient: Patient = None
        self._viewController = viewController

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self._dose1Label = QLabel('Dose 1:')
        self.layout.addWidget(self._dose1Label)
        self._dose1ComboBox = QComboBox(self)
        self.layout.addWidget(self._dose1ComboBox)

        self._dose2Label = QLabel('Dose 2:')
        self.layout.addWidget(self._dose2Label)
        self._dose2ComboBox = QComboBox(self)
        self.layout.addWidget(self._dose2ComboBox)

        self._runButton = QPushButton('Update!')
        self._runButton.clicked.connect(self._run)
        self.layout.addWidget(self._runButton)

        self.layout.addStretch()

        self.setCurrentPatient(self._viewController.currentPatient)
        self._viewController.currentPatientChangedSignal.connect(self.setCurrentPatient)

    def _selectedDose1(self):
        return self._dose1ComboBox.currentData()

    def _selectedDose2(self):
        return self._dose2ComboBox.currentData()

    def _fill_combo(self, combo: QComboBox, data_list, current_name=None):
        combo.clear()
        for d in data_list:
            combo.addItem(d.name, d)
        if current_name and data_list:
            names = [d.name for d in data_list]
            if current_name in names:
                combo.setCurrentIndex(names.index(current_name))
            else:
                combo.setCurrentIndex(0)

    def _update_dose_combos(self, _data=None):
        """Refresh Dose 1/Dose 2 lists. _data is the added/removed item (from signal)."""
        if self._patient is None:
            self._dose1ComboBox.clear()
            self._dose2ComboBox.clear()
            return
        data_list = _dose_like_data(self._patient)
        n1 = self._dose1ComboBox.currentText()
        n2 = self._dose2ComboBox.currentText()
        self._fill_combo(self._dose1ComboBox, data_list, n1)
        self._fill_combo(self._dose2ComboBox, data_list, n2)

    def setCurrentPatient(self, patient: Patient):
        if self._patient is not None:
            self._patient.patientDataAddedSignal.disconnect(self._update_dose_combos)
            self._patient.patientDataRemovedSignal.disconnect(self._update_dose_combos)
        self._patient = patient
        if self._patient is not None:
            self._patient.patientDataAddedSignal.connect(self._update_dose_combos)
            self._patient.patientDataRemovedSignal.connect(self._update_dose_combos)
        self._update_dose_combos()

    def _run(self):
        d1 = _to_dose_image(self._selectedDose1())
        d2 = _to_dose_image(self._selectedDose2())
        self._viewController.dose1 = d1
        self._viewController.dose2 = d2
