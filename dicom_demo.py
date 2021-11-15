import SimpleITK as sitk
import radiomics
from radiomics import featureextractor

reader = sitk.ImageSeriesReader()
dicom_names = reader.GetGDCMSeriesFileNames(path)
reader.SetFileNames(dicom_names)
image = reader.Execute()

