# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 13:41:08 2018

@author: 212710461
"""

import os
from radiomics import glcm, glrlm, imageoperations
import SimpleITK as sitk
import numpy as np
import scipy.misc
import nibabel as nib
import pydicom
import seaborn as sns
import warnings
from math import exp
from math import log
import logging
from pandas.api.types import is_string_dtype

# ---------------- ROI transformation --------------------------------

def translateROI(roi, step_size, tag):
    try:
        roi_array = sitk.GetArrayFromImage(roi)
        roi_array_temp = np.pad(roi_array,
                                ((0, 0),
                                 (abs(step_size), abs(step_size)),
                                 (abs(step_size), abs(step_size))),
                                'constant')
        if tag == 'left':
            roi_array_temp = roi_array_temp[:, abs(step_size):roi_array.shape[1] + abs(step_size), 2 * abs(step_size):]
        elif tag == 'right':
            roi_array_temp = roi_array_temp[:, abs(step_size):roi_array.shape[1] + abs(step_size), :roi_array.shape[2]]
        if tag == 'up':
            roi_array_temp = roi_array_temp[:, 2 * abs(step_size):, abs(step_size):roi_array.shape[2] + abs(step_size)]
        elif tag == 'down':
            roi_array_temp = roi_array_temp[:, :roi_array.shape[1], abs(step_size):roi_array.shape[2] + abs(step_size)]

        roi_translated = sitk.GetImageFromArray(roi_array_temp)
        roi_translated.CopyInformation(roi)
    except Exception as e:
        print(e)
        roi_translated = roi
    return roi_translated

def scale_roi(roi, meta_path, pixel_num, tag):
    '''

    :param roi: sitk.Image
    :param meta_path: String
    :param pixel_num: int
    :param tag: String
    :return: sitk.Image
    '''

    # Calculate how many minimeter needs to extend
    info_reader = sitk.ImageFileReader()
    info_reader.SetFileName(meta_path)
    info_reader.LoadPrivateTagsOn()
    info_reader.ReadImageInformation()
    mm_extended = pixel_num
    try:
        pixel_spacing = info_reader.GetMetaData('0028|0030')
        last_idx = pixel_spacing.rfind('\\')
        mm_extended = round(pixel_num / float(pixel_spacing[0:last_idx]))
    except Exception as e:
        logging.info("Scale ROI: " + e)
        logging.info('No pixel size exists. Using 1 mm/pixel.')

    try:
        roi_array = sitk.GetArrayFromImage(roi)
        roi_scaled = sitk.GetImageFromArray(roi_array)
        roi_scaled.CopyInformation(roi)
        for dummy in range(mm_extended):
            print("* ", end=' ')
            roi_array = sitk.GetArrayFromImage(roi_scaled)
            gif = sitk.GradientImageFilter()
            gradient_img = gif.Execute(roi_scaled)
            gradient_array = sitk.GetArrayFromImage(gradient_img)
            gradient_x = gradient_array[:, :, :, 0]
            gradient_y = gradient_array[:, :, :, 1]
            index_sets_x = np.where(gradient_x != 0)
            index_sets_y = np.where(gradient_y != 0)
            if tag == 'magnify':
                roi_array[index_sets_x] = 1
                roi_array[index_sets_y] = 1
            else:
                roi_array[index_sets_x] = 0
                roi_array[index_sets_y] = 0
            roi_scaled = sitk.GetImageFromArray(roi_array)
        roi_scaled.CopyInformation(roi)
    except Exception as e:
        logging.info(e)
    logging.info('Complete!')
    return (roi_scaled, None)

def scale_roi_and_save(tumor_image_roi_paths, tumor_image_roi, PaientID_list, pixel_num, tag, progressBar):

    peri_image_roi = {}
    peri_image_roi_paths = {}
    # extended_image_roi = {}
    # extended_image_roi_paths = {}

    lens = len(PaientID_list)
    step = 0

    for PatientID in PaientID_list:
        print('Scaling ' + PatientID + ' ...')
        image, roi, meta_path = tumor_image_roi[PatientID]
        image_path, roi_path, _ = tumor_image_roi_paths[PatientID]

        roi_scaled, _ = scale_roi(roi, meta_path, pixel_num, tag)

        # # Save extended ROI
        # path_scaled = './Extended_ROI/' + PatientID + '_ExtendedROI.nii'
        # writer = sitk.ImageFileWriter()
        # writer.SetFileName(path_scaled)
        # writer.Execute(roi_scaled)
        # img_original = nib.load(roi_path)
        # img_scaled = nib.load(path_scaled)
        # array_data = img_scaled.get_fdata()
        # array_data = array_data.astype(np.int16)
        # affine = img_original.affine
        # img3 = nib.Nifti1Image(array_data, affine)
        # img3.set_data_dtype(np.uint8)
        # nib.save(img3, path_scaled)
        # extended_roi = sitk.ReadImage(path_scaled)

        # Save peritumoral ROI
        roi_transformed_array = sitk.GetArrayFromImage(roi_scaled)
        roi_array = sitk.GetArrayFromImage(roi)
        roi_peri_array = np.absolute(roi_transformed_array - roi_array)
        roi_peri = sitk.GetImageFromArray(roi_peri_array)
        roi_peri.CopyInformation(roi)
        path_peri = './Peritumoral_ROI/' + PatientID + '_PeritumoralROI.nii'
        writer = sitk.ImageFileWriter()
        writer.SetFileName(path_peri)
        writer.Execute(roi_peri)

        img_original = nib.load(roi_path)
        img_peri = nib.load(path_peri)
        array_data = img_peri.get_fdata()
        array_data = array_data.astype(np.int16)
        affine = img_original.affine
        img3 = nib.Nifti1Image(array_data, affine)
        img3.set_data_dtype(np.uint8)
        nib.save(img3, path_peri)
        peri_roi = sitk.ReadImage(path_peri)

        peri_image_roi[PatientID] = (image, peri_roi, meta_path)
        peri_image_roi_paths[PatientID] = (image_path, path_peri, meta_path)
        # extended_image_roi[PatientID] = (image, extended_roi, meta_path)
        # extended_image_roi_paths[PatientID] = (image_path, path_scaled, meta_path)

        step += 1
        percentage = step / lens * 100
        progressBar.setValue(percentage)

    return (peri_image_roi, peri_image_roi_paths)
    # return (peri_image_roi, peri_image_roi_paths, extended_image_roi, extended_image_roi_paths)

# ---------------------------------------------------------------------



# ---------------Group Dicom folders and Nifti files -----------------------
def read_dicom_folder(path):
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(path)
    reader.SetFileNames(dicom_names)
    image = reader.Execute()
    return image

def grouping(path):
    dicom_files = {}
    nifti_files = {}
    for_meta = {}
    for root, dirs, files in os.walk(path):
        for name in files:
            each_dir = os.path.join(root, name)
            idx = each_dir.rfind('\\')
            # Grouping Dicom folders
            if not each_dir[0:idx] in dicom_files and each_dir.rfind('.nii') is -1:
                idx0 = each_dir.rfind('\\', 0, idx)
                name_node = each_dir[idx0 + 1:idx]
                dicom_files[name_node] = each_dir[0:idx]
                for_meta[name_node] = each_dir
            # Groupingf Nifti files
            else:
                idx1 = each_dir.rfind('_')
                name_node = each_dir[idx + 1:idx1]
                nifti_files[name_node] = each_dir

    return dicom_files, nifti_files, for_meta

def counts_to_suv(counts_image, meta_path):
    try:
        image_arr = sitk.GetArrayFromImage(counts_image)

        ds = pydicom.dcmread(meta_path)
        try:
            unit = ds[0x0054, 0x1001].value
        except:
            return ("Tag 0x0054, 0x1001 cannot be retrieved!")
        if unit == 'BQML':
            try:
                half_life = float(ds[0x0054, 0x0016][0][0x0018, 0x1075].value)
            except:
                return ("Tag 0x0018, 0x1075 cannot be retrieved!")
            try:
                scan_time = ds[0x0008, 0x0032].value[0:6]
            except:
                return ("Tag 0x0008, 0x0032 cannot be retrieved!")
            scan_time = int(scan_time[0:2]) * 3600 + \
                        int(scan_time[2:4]) * 60 + \
                        int(scan_time[4:6])
            try:
                start_time = ds[0x0054, 0x0016][0][0x0018, 0x1072].value[0:6]
            except:
                return ("Tag 0x0018, 0x1072 cannot be retrieved!")
            start_time = int(start_time[0:2]) * 3600 + \
                         int(start_time[2:4]) * 60 + \
                         int(start_time[4:6])
            decay_time = scan_time - start_time
            try:
                injected_dose = float(ds[0x0054, 0x0016][0][0x0018, 0x1074].value)
            except:
                return ("Tag 0x0018, 0x1074 cannot be retrieved!")
            decayed_dose = injected_dose * pow(2, -decay_time / half_life)
            try:
                weight = float(ds[0x0010, 0x1030].value)
            except:
                return ("Tag 0x0010, 0x1030 cannot be retrieved!")
            SUVbwScaleFactor = (weight * 1000 / decayed_dose)

        elif unit == 'CNTS':

            # Philips private scale factor
            try:
                SUVbwScaleFactor = float(ds[0x7053, 0x1000].value)
            except:
                try:
                    rescale_slope = float(ds[0x0028, 0x1053].value)
                except:
                    return ("Tag 0x0028, 0x1053 cannot be retrieved!")
                SUVbwScaleFactor = rescale_slope * float(ds[0x7053, 0x1009].value)

        elif unit == 'GML':

            SUVbwScaleFactor = 1.0
        try:
            rescale_slope = float(ds[0x0028, 0x1053].value)
        except:
            return ("Tag 0x0028, 0x1053 cannot be retrieved!")
        try:
            rescale_intercept = float(ds[0x0028, 0x1052].value)
        except:
            return ("Tag 0x0028, 0x1052 cannot be retrieved!")
        SUV_arr = (image_arr * rescale_slope + rescale_intercept) * SUVbwScaleFactor
        suv_image = sitk.GetImageFromArray(SUV_arr)
        suv_image.CopyInformation(counts_image)

        return suv_image

    except Exception as e:

        return e

def load_dicom_nii(folder_name, progressBar, suv_flag):
    '''
    Load dicom and nifti data in batch.
    Save the loaded images and ROIs in container.
    Save the corresponding paths of image folders and ROIs in container.
    Save the path of one of the dicom image for extracting header information in container.

    :return: [1] images_rois: list, each element contains (a, b, c, d), where a is image (sitk.Image),
                  b is ROI (sitk.Image), c is path of a dicom file for extracting meta data (String),
                  d is patient ID (String)
              [2] image_roi_paths: list, each element contain (a, b, c, d), where a is image folder path (String),
                  b is ROI file path (String), c is path of a dicom file for extracting meta data(String),
                  d is patient ID (String)
              [3] PatientID_list: list, store a list of Patient ID
    '''
    error_roi_list = []
    no_matching_nifti = []
    PatientID_list = []
    images_rois = {}
    image_roi_paths = {}
    dicom_files, nifti_files, for_meta = grouping(folder_name)
    lens = len(dicom_files)
    step = 0
    for PatientID in dicom_files.keys():
        if PatientID in nifti_files.keys():

            try:
                logging.info('Loading ' + PatientID + ' ...')

                # Extract dicom's and nifti's path
                image_path = dicom_files[PatientID]
                roi_path = nifti_files[PatientID]
                meta_path = for_meta[PatientID]

                image = read_dicom_folder(image_path)
                roi = sitk.ReadImage(roi_path)

                # Check whether ROI is OK to be used
                label_arr = sitk.GetArrayFromImage(roi)
                label_info = np.unique(label_arr)
                if len(label_info) != 2:
                    error_roi_list.append(PatientID + '    -> More than two values!')
                elif 0 not in label_info or 1 not in label_info:
                    # print(label_info, label_info[1], type(label_info[1]))
                    num = str(label_info[1])
                    logging.info('The original label value is ' + num +  ', change the label value into 1')
                    label_arr -= (label_info[1] - 1)
                    roi = sitk.GetImageFromArray(label_arr)
                    roi.CopyInformation(sitk.ReadImage(roi_path))
                    # error_roi_list.append(PatientID + '    -> The label values do not contain 0 or 1')

                # Check if PET image, yes then transform SUV
                if suv_flag:
                    info_reader = sitk.ImageFileReader()
                    info_reader.SetFileName(meta_path)
                    info_reader.LoadPrivateTagsOn()
                    info_reader.ReadImageInformation()
                    try:
                        modality = info_reader.GetMetaData('0008|0060')
                        logging.info(modality)
                        corrected_image = info_reader.GetMetaData('0028|0051')
                        logging.info(corrected_image)
                        decay = info_reader.GetMetaData('0054|1102')
                        logging.info(decay)

                        if modality == 'PT' and 'ATTN' in corrected_image and 'DECY' in corrected_image and 'START' in decay:
                            image = counts_to_suv(image, meta_path)

                        if type(image) == str:
                            logging.error(image)
                            logging.error(PatientID + ' needs to be removed!')
                            continue
                        else:
                            logging.info('Counts image transformed to SUV image!')
                    except:
                        modality = None
                        logging.info('No enough info for converting SUV ', decay, corrected_image)

                # Check if RGB then change to gray scale
                image_array = sitk.GetArrayFromImage(image)
                if len(image_array.shape) > 3:
                    gray_image_array = sitk.GetArrayFromImage(image)[:, :, :, 0] * 0.2126 + \
                                       sitk.GetArrayFromImage(image)[:, :, :, 1] * 0.7152 + \
                                       sitk.GetArrayFromImage(image)[:, :, :, 2] * 0.0722

                    gray_image_array = gray_image_array.astype('uint8')
                    gray_image = sitk.GetImageFromArray(gray_image_array)
                    gray_image.CopyInformation(roi)
                    image = gray_image

                settings = {'correctMask': True,
                            'geometryTolerance': 100}


                # Check whether ROI is marching
                if imageoperations.checkMask(image, roi, **settings)[1]:
                    logging.info('ROI corrected!')
                    roi = imageoperations.checkMask(image, roi, **settings)[1]

                PatientID_list.append(PatientID)
                images_rois[PatientID] = (image, roi, meta_path)
                image_roi_paths[PatientID] = (image_path, roi_path, meta_path)

                logging.info('Done!')

            except Exception as e:
                logging.info('Load ' + PatientID + ' failed!')
                logging.info(e)
                no_matching_nifti.append(PatientID)
        else:
            no_matching_nifti.append(PatientID)

        step += 1
        percentage = step / lens * 100
        progressBar.setValue(percentage)

    if no_matching_nifti != []:
        logging.info('The following nifti files failed the matching. ' + '\n' + str(no_matching_nifti))

    return (images_rois, image_roi_paths, PatientID_list, error_roi_list)

# ----------------------------------------------------------------------------


# ------------------------------- Obtain images and plot them ---------------------------------
def plot_image_info(subplot, image, mask, gs):

    warnings.filterwarnings('ignore')

    # Plot image with ROI
    img = sitk.Cast(sitk.RescaleIntensity(image), sitk.sitkUInt8)
    # seg = sitk.Cast(sitk.RescaleIntensity(mask), sitk.sitkUInt8)
    seg = mask
    img_arr = sitk.GetArrayFromImage(img)
    roi_arr = sitk.GetArrayFromImage(seg)

    sum = []
    for i in range(roi_arr.shape[0]):
        sum.append(roi_arr[i].sum())
    idx = sum.index(max(sum))

    img = sitk.Image(3, img_arr.shape[2], img_arr.shape[1], sitk.sitkUInt8)
    overlay = sitk.GetArrayFromImage(img)
    overlay[:, :, 0] = img_arr[idx]
    overlay[:, :, 1] = img_arr[idx] - img_arr[idx] * roi_arr[idx] * 0.5
    overlay[:, :, 2] = img_arr[idx] - img_arr[idx] * roi_arr[idx] * 0.5

    ax = subplot.add_subplot(gs[:, :3])
    ax.axis('off')
    ax.set_title('ROI')
    ax.imshow(overlay)

    # Plot the GLCM
    cal_glcm = glcm.RadiomicsGLCM(image, mask)
    cal_glcm.settings['distance'] = [1]
    glcm_all_angle = cal_glcm._calculateMatrix()
    glcm_mean = np.mean(glcm_all_angle, axis=2)
    glcm_mean = scipy.misc.imresize(glcm_mean, size=(25, 25))
    ax = subplot.add_subplot(gs[0, 3])
    ax.set_title('GLCM', fontsize='x-small')
    sns.heatmap(glcm_mean, ax=ax)
    ax.tick_params(bottom=False)
    ax.tick_params(labelsize='x-small')

    # Plot GLRLM
    cal_rlm = glrlm.RadiomicsGLRLM(image, mask)
    cal_rlm._initCalculation()
    rlm_mean = cal_rlm._calculateMatrix()
    rlm_mean = np.mean(rlm_mean, axis=2)[:, 0:25]
    rlm_mean = scipy.misc.imresize(rlm_mean, size=(25, 25))
    ax = subplot.add_subplot(gs[1, 3])
    ax.set_title('RLM', fontsize='x-small')
    sns.heatmap(rlm_mean, ax=ax)
    ax.tick_params(bottom=False)
    ax.tick_params(labelsize='x-small')

    # Plot histogram
    bbox, _ = imageoperations.checkMask(image, mask)
    cropImage, _ = imageoperations.cropToTumorMask(image, mask, bbox)
    cropImage_arr = np.squeeze(sitk.GetArrayFromImage(cropImage))
    cropImage_arr = cropImage_arr.flatten()
    ax = subplot.add_subplot(gs[2, 3])
    ax.autoscale()
    ax.hist(cropImage_arr, bins=256)
    ax.set_title('Histogram')
    ax.tick_params(length=0.2)

    subplot.set_size_inches(8, 6)
    subplot.tight_layout()

    warnings.filterwarnings('default')

    return (0, 0)

def plot_scaled_roi(subplot, image, mask):

    try:
        # Plot image with ROI
        img = sitk.Cast(sitk.RescaleIntensity(image), sitk.sitkUInt8)
        # seg = sitk.Cast(sitk.RescaleIntensity(mask), sitk.sitkUInt8)
        seg = mask
        img_arr = sitk.GetArrayFromImage(img)
        roi_arr = sitk.GetArrayFromImage(seg)

        sum = []
        for i in range(roi_arr.shape[0]):
            sum.append(roi_arr[i].sum())
        idx = sum.index(max(sum))

        img = sitk.Image(3, img_arr.shape[2], img_arr.shape[1], sitk.sitkUInt8)
        overlay = sitk.GetArrayFromImage(img)
        overlay[:, :, 0] = img_arr[idx]
        overlay[:, :, 1] = img_arr[idx] - img_arr[idx] * roi_arr[idx] * 0.5
        overlay[:, :, 2] = img_arr[idx] - img_arr[idx] * roi_arr[idx] * 0.5

        ax = subplot.add_subplot(111)
        ax.axis('off')
        ax.set_title('Transformed ROI')
        ax.imshow(overlay)

        subplot.tight_layout()

    except Exception as e:
        logging.info(e)

    return (0, 0)
# --------------------------------------------------------------------------------------


# ============ Read image and roi ======================================
def readImageROI(image_path, roi_path, meta_path):

    image = read_dicom_folder(image_path)
    label = sitk.ReadImage(roi_path)

    # Check whether ROI is OK to be used
    label_arr = sitk.GetArrayFromImage(label)
    label_info = np.unique(label_arr)
    if 0 not in label_info or 1 not in label_info:
        # print(label_info, label_info[1], type(label_info[1]))
        num = str(label_info[1])
        logging.info('The original label value is ' + num + ', change the label value into 1')
        label_arr -= (label_info[1] - 1)
        label = sitk.GetImageFromArray(label_arr)
        label.CopyInformation(sitk.ReadImage(roi_path))
        # error_roi_list.append(PatientID + '    -> The label values do not contain 0 or 1')

    # Check if PET image, yes then transform SUV
    info_reader = sitk.ImageFileReader()
    info_reader.SetFileName(meta_path)
    info_reader.LoadPrivateTagsOn()
    info_reader.ReadImageInformation()
    try:
        modality = info_reader.GetMetaData('0008|0060')
        corrected_image = info_reader.GetMetaData('0028|0051')
        decay = info_reader.GetMetaData('0054|1102')
        if modality == 'PT' and 'ATTN' in corrected_image and 'DECY' in corrected_image and 'START' in decay:
            image = counts_to_suv(image, meta_path)
    except:
        modality = None


    # Check if RGB then change to gray scale
    image_array = sitk.GetArrayFromImage(image)
    if len(image_array.shape) > 3:
        gray_image_array = sitk.GetArrayFromImage(image)[:, :, :, 0] * 0.2126 + \
                           sitk.GetArrayFromImage(image)[:, :, :, 1] * 0.7152 + \
                           sitk.GetArrayFromImage(image)[:, :, :, 2] * 0.0722

        gray_image_array = gray_image_array.astype('uint8')
        gray_image = sitk.GetImageFromArray(gray_image_array)
        gray_image.CopyInformation(label)
        image = gray_image

    settings = {'correctMask': True,
                'geometryTolerance': 100}

    if imageoperations.checkMask(image, label, **settings)[1]:
        logging.info('Mask corrected!')
        label = imageoperations.checkMask(image, label, **settings)[1]

    boundingBox = imageoperations.checkMask(image, label, **settings)[0]

    extracted_image, extracted_label = imageoperations.cropToTumorMask(image, label, boundingBox)

    return extracted_image, extracted_label

def readImageROI_noCrop(image_path, roi_path, meta_path):

    image = read_dicom_folder(image_path)
    label = sitk.ReadImage(roi_path)

    # Check whether ROI is OK to be used
    label_arr = sitk.GetArrayFromImage(label)
    label_info = np.unique(label_arr)
    if 0 not in label_info or 1 not in label_info:
        # print(label_info, label_info[1], type(label_info[1]))
        num = str(label_info[1])
        logging.info('The original label value is ' + num + ', change the label value into 1')
        label_arr -= (label_info[1] - 1)
        label = sitk.GetImageFromArray(label_arr)
        label.CopyInformation(sitk.ReadImage(roi_path))
        # error_roi_list.append(PatientID + '    -> The label values do not contain 0 or 1')

    # Check if PET image, yes then transform SUV
    info_reader = sitk.ImageFileReader()
    info_reader.SetFileName(meta_path)
    info_reader.LoadPrivateTagsOn()
    info_reader.ReadImageInformation()
    try:
        modality = info_reader.GetMetaData('0008|0060')
        corrected_image = info_reader.GetMetaData('0028|0051')
        decay = info_reader.GetMetaData('0054|1102')

        if modality == 'PT' and 'ATTN' in corrected_image and 'DECY' in corrected_image and 'START' in decay:
            image = counts_to_suv(image, meta_path)
    except:
        modality = None

    # Check if RGB then change to gray scale
    image_array = sitk.GetArrayFromImage(image)
    if len(image_array.shape) > 3:
        gray_image_array = sitk.GetArrayFromImage(image)[:, :, :, 0] * 0.2126 + \
                           sitk.GetArrayFromImage(image)[:, :, :, 1] * 0.7152 + \
                           sitk.GetArrayFromImage(image)[:, :, :, 2] * 0.0722

        gray_image_array = gray_image_array.astype('uint8')
        gray_image = sitk.GetImageFromArray(gray_image_array)
        gray_image.CopyInformation(label)
        image = gray_image

    settings = {'correctMask': True,
                'geometryTolerance': 100}

    if imageoperations.checkMask(image, label, **settings)[1]:
        logging.info('Mask corrected!')
        label = imageoperations.checkMask(image, label, **settings)[1]

    return image, label

# ======================================================================

def check_string(df):

    string_list = []
    columns_list = df.columns.tolist()
    for column in columns_list:
        if is_string_dtype(df[column]):
            string_list.append(column)

    return (string_list, df)