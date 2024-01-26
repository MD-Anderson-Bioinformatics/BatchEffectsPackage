#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright (c) 2011-2024 University of Texas MD Anderson Cancer Center

This program is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 2 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.

MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>
@author: Tod Casasent
"""


import io
import os
import zipfile
import json
from typing import Dict, List
import pandas
import xmltodict
from mbatch.test.common import xml_text, write_tsv_list, get_entry_as_list
from mbatch.test.common import get_data_string_text
from mbatch.test.common import get_tsv_entry_dict, get_dataframe_entry
from mbatch.gdcapi.download_clinical import GdcApiClinical
from mbatch.test.common import add_error

CLINICAL_HEADERS: List[str] = [
    'bcr_patient_barcode',
    'bcr_patient_uuid',
    'tumor_site',
    'histologic_diagnosis',
    'history_other_malignancy',
    'sex',
    'vital_status',
    'days_to_birth',
    # 'days_to_last_known_alive',
    'days_to_death',
    'days_to_last_followup',
    'history_of_neoadjuvant_treatment',
    'ethnicity',
    'height',
    'weight',
    'person_neoplasm_cancer_status',
    'days_to_initial_pathologic_diagnosis',
    'age_at_initial_pathologic_diagnosis',
    'year_of_initial_pathologic_diagnosis',
    'tumor_type',
    'tobacco_smoking_history_indicator',
    'tobacco_smoking_year_started',
    'tobacco_smoking_year_stopped',
    'tobacco_smoking_pack_years_smoked',
    'pharmaceutical_tx_adjuvant',
    'radiation_treatment_adjuvant',
    'treatment_outcome_first_course',
    'race'
]


# pylint: disable=too-many-instance-attributes,too-many-arguments,too-few-public-methods
class Clinical:
    """
    Class to hold public clinical data for a patient.
    """
    # declare but do not set member attributes
    bcr_patient_barcode: str
    bcr_patient_uuid: str
    tumor_site: str
    histologic_diagnosis: str
    history_other_malignancy: str
    sex: str
    vital_status: str
    days_to_birth: str
    days_to_last_known_alive: str
    days_to_death: str
    days_to_last_followup: str
    history_of_neoadjuvant_treatment: str
    ethnicity: str
    height: str
    weight: str
    person_neoplasm_cancer_status: str
    days_to_initial_pathologic_diagnosis: str
    age_at_initial_pathologic_diagnosis: str
    year_of_initial_pathologic_diagnosis: str
    tumor_type: str
    tobacco_smoking_history_indicator: str
    tobacco_smoking_year_started: str
    tobacco_smoking_year_stopped: str
    tobacco_smoking_pack_years_smoked: str
    pharmaceutical_tx_adjuvant: str
    radiation_treatment_adjuvant: str
    treatment_outcome_first_course: str
    race: str

    def __init__(self: 'Clinical', the_dict: Dict[str, str]) -> None:
        """
        Copy dictionary values into attributes.
        :param the_dict: Dictionary where keys match attribute names of Clinical class.
        """
        super().__init__()
        # declare instance variables
        self.bcr_patient_barcode: str = 'Unknown'
        self.bcr_patient_uuid: str = 'Unknown'
        self.tumor_site: str = 'Unknown'
        self.histologic_diagnosis: str = 'Unknown'
        self.history_other_malignancy: str = 'Unknown'
        self.sex: str = 'Unknown'
        self.vital_status: str = 'Unknown'
        self.days_to_birth: str = 'Unknown'
        self.days_to_last_known_alive: str = 'Unknown'
        self.days_to_death: str = 'Unknown'
        self.days_to_last_followup: str = 'Unknown'
        self.history_of_neoadjuvant_treatment: str = 'Unknown'
        self.ethnicity: str = 'Unknown'
        self.height: str = 'Unknown'
        self.weight: str = 'Unknown'
        self.person_neoplasm_cancer_status: str = 'Unknown'
        self.days_to_initial_pathologic_diagnosis: str = 'Unknown'
        self.age_at_initial_pathologic_diagnosis: str = 'Unknown'
        self.year_of_initial_pathologic_diagnosis: str = 'Unknown'
        self.tumor_type: str = 'Unknown'
        self.tobacco_smoking_history_indicator: str = 'Unknown'
        self.tobacco_smoking_year_started: str = 'Unknown'
        self.tobacco_smoking_year_stopped: str = 'Unknown'
        self.tobacco_smoking_pack_years_smoked: str = 'Unknown'
        self.pharmaceutical_tx_adjuvant: str = 'Unknown'
        self.radiation_treatment_adjuvant: str = 'Unknown'
        self.treatment_outcome_first_course: str = 'Unknown'
        self.race: str = 'Unknown'
        # set values
        attr: str
        for attr in CLINICAL_HEADERS:
            setattr(self, attr, the_dict[attr])

    def write_clinical(self: 'Clinical', the_out_file: io.TextIOWrapper) -> None:
        """
        Write out just the clinical portion of a file element to a file,
        after first calling parent history to write out its portion.
        Writes a newline.
        :param the_out_file: stream to which to write strings
        :return: None
        """
        data_list: List[str] = []
        for attr in CLINICAL_HEADERS:
            data_list.append(getattr(self, attr))
        write_tsv_list(the_out_file, data_list, True, False)
# pylint: enable=too-many-instance-attributes,too-many-arguments,too-few-public-methods

# ########################################################
# writing files
# ########################################################


def add_clinical_file_clinical_xml(the_archive_file: str, the_xml_file: str,
                                   the_clinical_info: Dict[str, Clinical]) -> None:
    """
    Read Public Clinical XML file and convert to Dictionary of Clinical objects.
    :param the_archive_file: the full path and filename for archive (ZIP) file
    :param the_xml_file: the name of the file inside the ZIP file to actually read
    :param the_clinical_info: Dictionary of Clinical objects (key may be UUID or Barcode, depending on creation)
    :return: Nothing
    """
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    with zipfile.ZipFile(the_archive_file, 'r') as zip_file:
        with zip_file.open(the_xml_file, mode="r") as in_file:
            # need type
            doc: Dict = xmltodict.parse(in_file.read())
            if 'ssf:tcga_bcr' in doc:
                print("add_clinical_file skip ssf file", flush=True)
            else:
                # get base key, which should be XXXX:tcga_bcr
                base_key: str = list(doc.keys())[0]
                if not base_key.endswith(":tcga_bcr"):
                    raise ValueError(f"Expected XXXX:tcga_bcr but found {base_key}")
                base_prefix: str = base_key.removesuffix(":tcga_bcr")
                patient: Dict = doc[base_key][base_prefix+':patient']
                clinical_dict: Dict[str, str] = \
                    {
                        'bcr_patient_barcode': get_data_string_text(patient, 'shared:bcr_patient_barcode', 'Unknown'),
                        'bcr_patient_uuid': get_data_string_text(patient, 'shared:bcr_patient_uuid', 'Unknown'),
                        'tumor_site': get_data_string_text(patient, 'clin_shared:tumor_tissue_site', 'Unknown'),
                        'histologic_diagnosis': get_data_string_text(patient, 'shared:histological_type', 'Unknown'),
                        'history_other_malignancy': get_data_string_text(patient, 'shared:other_dx', 'Unknown'),
                        'sex': get_data_string_text(patient, 'shared:gender', 'Unknown'),
                        'vital_status': get_data_string_text(patient, 'clin_shared:vital_status', 'Unknown'),
                        'days_to_birth': get_data_string_text(patient, 'clin_shared:days_to_birth', 'Unknown'),
                        # 'days_to_last_known_alive': get_data_string_text(patient, 'clin_shared:days_to_last_known_alive', 'Unknown'),
                        'days_to_death': get_data_string_text(patient, 'clin_shared:days_to_death', 'Unknown'),
                        'days_to_last_followup': get_data_string_text(patient, 'clin_shared:days_to_last_followup', 'Unknown'),
                        'history_of_neoadjuvant_treatment': get_data_string_text(patient, 'shared:history_of_neoadjuvant_treatment', 'Unknown'),
                        'ethnicity': get_data_string_text(patient, 'clin_shared:ethnicity', 'Unknown'),
                        'height': get_data_string_text(patient, 'clin_shared:height', 'Unknown'),
                        'weight': get_data_string_text(patient, 'clin_shared:weight', 'Unknown'),
                        'person_neoplasm_cancer_status': get_data_string_text(patient, 'clin_shared:person_neoplasm_cancer_status', 'Unknown'),
                        'days_to_initial_pathologic_diagnosis': get_data_string_text(patient, 'clin_shared:days_to_initial_pathologic_diagnosis', 'Unknown'),
                        'age_at_initial_pathologic_diagnosis': get_data_string_text(patient, 'clin_shared:age_at_initial_pathologic_diagnosis', 'Unknown'),
                        'year_of_initial_pathologic_diagnosis': get_data_string_text(patient, 'clin_shared:year_of_initial_pathologic_diagnosis', 'Unknown'),
                        'tumor_type': get_data_string_text(patient, base_prefix + ':tumor_type', 'Unknown'),
                        'tobacco_smoking_history_indicator': get_data_string_text(patient, 'shared:tobacco_smoking_history', 'Unknown'),
                        'tobacco_smoking_year_started': get_data_string_text(patient, 'clin_shared:year_of_tobacco_smoking_onset', 'Unknown'),
                        'tobacco_smoking_year_stopped': get_data_string_text(patient, 'clin_shared:stopped_smoking_year', 'Unknown'),
                        'tobacco_smoking_pack_years_smoked': get_data_string_text(patient, 'clin_shared:number_pack_years_smoked', 'Unknown'),
                        'pharmaceutical_tx_adjuvant': get_data_string_text(patient, 'clin_shared:postoperative_rx_tx', 'Unknown'),
                        'radiation_treatment_adjuvant': get_data_string_text(patient, 'clin_shared:radiation_therapy', 'Unknown'),
                        'treatment_outcome_first_course': get_data_string_text(patient, 'clin_shared:primary_therapy_outcome_success', 'Unknown')
                    }
                #
                race_list: List[Dict] = get_entry_as_list(patient, 'clin_shared:race_list', 'clin_shared:race')
                race: str = ''
                val: Dict
                for val in race_list:
                    if '' == race:
                        race = xml_text(val, 'Unknown')
                    else:
                        race = race + "," + xml_text(val, 'Unknown')
                clinical_dict['race'] = race
                b_clinical: Clinical = Clinical(clinical_dict)
                the_clinical_info[b_clinical.bcr_patient_uuid] = b_clinical


def add_clinical_file_clinical_tsv(the_archive_file: str, the_tsv_file: str,
                                   the_clinical_info: Dict[str, Clinical]) -> None:
    """
    Read Public Clinical TSV file and convert to Dictionary of Clinical objects.
    :param the_archive_file: the full path and filename for archive (ZIP) file
    :param the_tsv_file: the name of the file inside the ZIP file to actually read
    :param the_clinical_info: Dictionary of Clinical objects (key may be UUID or Barcode, depending on creation)
    :return: Nothing
    """
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    with zipfile.ZipFile(the_archive_file, 'r') as zip_file:
        with zip_file.open(the_tsv_file, mode="r") as in_file:
            index: int = 0
            keys: List[str]
            for line in io.TextIOWrapper(in_file, encoding="utf-8"):
                if 0 == index % 1000:
                    print(f"add_clinical_file_clinical_tsv datafile index={index}", flush=True)
                if 0 == index:
                    keys = line.rstrip('\n').split('\t')
                else:
                    values: List[str] = line.rstrip('\n').split("\t")
                    patient_line: Dict[str, str] = dict(zip(keys, values))
                    clinical_dict: Dict[str, str] = \
                        {
                            'bcr_patient_barcode': get_tsv_entry_dict(patient_line, 'case_id', 'Unknown'),
                            'bcr_patient_uuid': get_tsv_entry_dict(patient_line, 'cases.submitter_id', 'Unknown'),
                            'tumor_site': get_tsv_entry_dict(patient_line, 'diagnoses.site_of_resection_or_biopsy', 'Unknown'),
                            'histologic_diagnosis': get_tsv_entry_dict(patient_line, 'diagnoses.classification_of_tumor', 'Unknown'),
                            'history_other_malignancy': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'sex': get_tsv_entry_dict(patient_line, 'demographic.gender', 'Unknown'),
                            'race': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'vital_status': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'days_to_birth': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'days_to_death': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'days_to_last_followup': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'history_of_neoadjuvant_treatment': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'ethnicity': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'height': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'weight': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'person_neoplasm_cancer_status': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'days_to_initial_pathologic_diagnosis': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'age_at_initial_pathologic_diagnosis': get_tsv_entry_dict(patient_line, 'diagnoses.age_at_diagnosis', 'Unknown'),
                            'year_of_initial_pathologic_diagnosis': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'tumor_type': get_tsv_entry_dict(patient_line, 'diagnoses.morphology', 'Unknown'),
                            'tobacco_smoking_history_indicator': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'tobacco_smoking_year_started': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'tobacco_smoking_year_stopped': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'tobacco_smoking_pack_years_smoked': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'pharmaceutical_tx_adjuvant': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'radiation_treatment_adjuvant': get_tsv_entry_dict(patient_line, '', 'Unknown'),
                            'treatment_outcome_first_course': get_tsv_entry_dict(patient_line, '', 'Unknown')
                        }
                    b_clinical: Clinical = Clinical(clinical_dict)
                    the_clinical_info[b_clinical.bcr_patient_uuid] = b_clinical
                index += 1


# pylint: disable=too-many-nested-blocks,too-many-branches,too-many-locals,too-many-statements
def add_clinical_file_clinical_json(the_archive_file: str, the_json_file: str,
                                    the_clinical_info: Dict[str, Clinical]) -> None:
    """
    Read Public Clinical JSON file and convert to Dictionary of Clinical objects.
    :param the_archive_file: the full path and filename for archive (ZIP) file
    :param the_json_file: the name of the file inside the ZIP file to actually read
    :param the_clinical_info: Dictionary of Clinical objects (key may be UUID or Barcode, depending on creation)
    :return: Nothing
    """
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    with zipfile.ZipFile(the_archive_file, 'r') as zip_file:
        with zip_file.open(the_json_file, mode="r") as in_file:
            json_obj: Dict = json.loads(in_file.read())
            clinical: Dict[str, str] = \
                {
                    'bcr_patient_barcode': json_obj['ClinicalData']['SubjectData']['submitter_id'],
                    'bcr_patient_uuid': json_obj['ClinicalData']['SubjectData']['id'],
                    'tumor_site': 'Unknown', 'histologic_diagnosis': 'Unknown', 'history_other_malignancy': 'Unknown',
                    'sex': 'Unknown', 'vital_status': 'Unknown', 'days_to_birth': 'Unknown', 'days_to_death': 'Unknown',
                    'days_to_last_followup': 'Unknown', 'history_of_neoadjuvant_treatment': 'Unknown',
                    'ethnicity': 'Unknown', 'height': 'Unknown', 'weight': 'Unknown',
                    'person_neoplasm_cancer_status': 'Unknown', 'days_to_initial_pathologic_diagnosis': 'Unknown',
                    'age_at_initial_pathologic_diagnosis': 'Unknown', 'year_of_initial_pathologic_diagnosis': 'Unknown',
                    'tumor_type': 'Unknown', 'tobacco_smoking_history_indicator': 'Unknown', 'tobacco_smoking_year_started': 'Unknown',
                    'tobacco_smoking_year_stopped': 'Unknown', 'tobacco_smoking_pack_years_smoked': 'Unknown',
                    'pharmaceutical_tx_adjuvant': 'Unknown', 'radiation_treatment_adjuvant': 'Unknown',
                    'treatment_outcome_first_course': 'Unknown', 'race': 'Unknown'
                }
            entity_list: List[Dict] = json_obj['ClinicalData']['GDC Clinical Data Entities']
            entity: Dict
            for entity in entity_list:
                entity_type: str = entity['type']
                if 'case' == entity_type:
                    clinical['tumor_site'] = entity.get('primary_site', 'Unknown')
                elif 'demographic' == entity_type:
                    clinical['sex'] = entity.get('gender', 'Unknown')
                    clinical['vital_status'] = entity.get('vital_status', 'Unknown')
                    clinical['days_to_birth'] = entity.get('days_to_birth', 'Unknown')
                    clinical['ethnicity'] = entity.get('ethnicity', 'Unknown')
                    clinical['race'] = entity.get('race', 'Unknown')
                elif 'diagnosis' == entity_type:
                    clinical['histologic_diagnosis'] = entity.get('ajcc_clinical_stage', 'Unknown')
                    clinical['tumor_type'] = entity.get('tumor_grade', 'Unknown')
                    clinical['history_other_malignancy'] = entity.get('prior_malignancy', 'Unknown')
                elif 'exposure' == entity_type:
                    clinical['tobacco_smoking_history_indicator'] = entity.get('tobacco_smoking_status', 'Unknown')
                elif 'follow_up' == entity_type:
                    value: str = entity['days_to_follow_up']
                    if value is not None:
                        if 'Unknown' == clinical['days_to_last_followup']:
                            clinical['days_to_last_followup'] = value
                        elif float(value) > float(clinical['days_to_last_followup']):
                            clinical['days_to_last_followup'] = value
                        else:
                            clinical['days_to_last_followup'] = value
                elif ('treatment' == entity_type) & ('Neoadjuvant' == entity.get('treatment_intent_type', 'Unknown')):
                    value: str = entity['treatment_or_therapy']
                    if 'no' == clinical['days_to_last_followup']:
                        clinical['history_of_neoadjuvant_treatment'] = value
                    else:
                        clinical['history_of_neoadjuvant_treatment'] = value
                elif ('treatment' == entity_type) & ('Adjuvant' == entity.get('treatment_intent_type', 'Unknown')) & \
                        ('Pharmaceutical' in entity.get('treatment_type', 'Unknown')):
                    value: str = entity['treatment_or_therapy']
                    if 'no' == clinical['days_to_last_followup']:
                        clinical['pharmaceutical_tx_adjuvant'] = value
                    else:
                        clinical['pharmaceutical_tx_adjuvant'] = value
                elif ('treatment' == entity_type) & ('Adjuvant' == entity.get('treatment_intent_type', 'Unknown')) & \
                        ('Radiation' in entity.get('treatment_type', 'Unknown')):
                    value: str = entity['treatment_or_therapy']
                    if 'no' == clinical['radiation_treatment_adjuvant']:
                        clinical['radiation_treatment_adjuvant'] = value
                    else:
                        clinical['radiation_treatment_adjuvant'] = value
            study_event_list: List[Dict] = json_obj['ClinicalData']['SubjectData']['StudyEventData']
            study_event: Dict
            for study_event in study_event_list:
                form_data_list: List[Dict] = study_event['FormData']
                form_data: Dict
                for form_data in form_data_list:
                    item_group_data_list: List[Dict] = form_data['ItemGroupData']
                    item_group: Dict
                    for item_group in item_group_data_list:
                        item_data_list: List[Dict] = item_group['ItemData']
                        item_data: Dict
                        for item_data in item_data_list:
                            if 'AgeAtDeath' == item_data['@ItemOID']:
                                clinical['days_to_death'] = item_data['@Value']
                            elif 'Height' == item_data['@ItemOID']:
                                clinical['height'] = item_data['@Value']
                            elif 'Weight' == item_data['@ItemOID']:
                                clinical['weight'] = item_data['@Value']
            my_clinical: Clinical = Clinical(clinical)
            the_clinical_info[my_clinical.bcr_patient_uuid] = my_clinical
# pylint: enable=too-many-nested-blocks,too-many-branches,too-many-locals,too-many-statements


def add_clinical_file_clinical_xlsx(the_archive_file: str, the_xlsx_file: str,
                                    the_clinical_info: Dict[str, Clinical]) -> None:
    """
    Read Public Clinical XLSX file and convert to Dictionary of Clinical objects.
    :param the_archive_file: the full path and filename for archive (ZIP) file
    :param the_xlsx_file: the name of the file inside the ZIP file to actually read
    :param the_clinical_info: Dictionary of Clinical objects (key may be UUID or Barcode, depending on creation)
    :return: Nothing
    """
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    with zipfile.ZipFile(the_archive_file) as zip_file:
        with zip_file.open(the_xlsx_file) as in_file:
            # add ,0 to target desired read_excel version
            my_df: pandas.DataFrame = pandas.read_excel(in_file, 0)
            # convert missing values to empty string
            my_df = my_df.mask(my_df.isna(), '')
            # iterate rows of dataframe
            index: int
            row: pandas.Series
            for index, row in my_df.iterrows():
                if 0 == index % 1000:
                    print(f"add_clinical_file_clinical_xlsx index={index}", flush=True)
                clinical_dict: Dict[str, str] = \
                    {
                        'bcr_patient_barcode': get_dataframe_entry(row, 'TARGET USI', 'Unknown'),
                        'bcr_patient_uuid': get_dataframe_entry(row, '', 'Unknown'),
                        'tumor_site': get_dataframe_entry(row, '', 'Unknown'),
                        'histologic_diagnosis': get_dataframe_entry(row, '', 'Unknown'),
                        'history_other_malignancy': get_dataframe_entry(row, '', 'Unknown'),
                        'sex': get_dataframe_entry(row, 'Gender', 'Unknown'),
                        'race': get_dataframe_entry(row, 'Race', 'Unknown'),
                        'vital_status': get_dataframe_entry(row, 'Vital Status', 'Unknown'),
                        'days_to_birth': get_dataframe_entry(row, 'Overall Survival Time in Days', 'Unknown'),
                        'days_to_death': get_dataframe_entry(row, '', 'Unknown'),
                        'days_to_last_followup': get_dataframe_entry(row, '', 'Unknown'),
                        'history_of_neoadjuvant_treatment': get_dataframe_entry(row, '', 'Unknown'),
                        'ethnicity': get_dataframe_entry(row, 'Ethnicity', 'Unknown'),
                        'height': get_dataframe_entry(row, '', 'Unknown'),
                        'weight': get_dataframe_entry(row, '', 'Unknown'),
                        'person_neoplasm_cancer_status': get_dataframe_entry(row, '', 'Unknown'),
                        'days_to_initial_pathologic_diagnosis': get_dataframe_entry(row, '', 'Unknown'),
                        'age_at_initial_pathologic_diagnosis': get_dataframe_entry(row, 'Age at Diagnosis in Days', 'Unknown'),
                        'year_of_initial_pathologic_diagnosis': get_dataframe_entry(row, '', 'Unknown'),
                        'tumor_type': get_dataframe_entry(row, 'FAB Category', 'Unknown'),
                        'tobacco_smoking_history_indicator': get_dataframe_entry(row, '', 'Unknown'),
                        'tobacco_smoking_year_started': get_dataframe_entry(row, '', 'Unknown'),
                        'tobacco_smoking_year_stopped': get_dataframe_entry(row, '', 'Unknown'),
                        'tobacco_smoking_pack_years_smoked': get_dataframe_entry(row, '', 'Unknown'),
                        'pharmaceutical_tx_adjuvant': get_dataframe_entry(row, 'Gemtuzumab ozogamicin treatment', 'Unknown'),
                        'radiation_treatment_adjuvant': get_dataframe_entry(row, '', 'Unknown'),
                        'treatment_outcome_first_course': get_dataframe_entry(row, '', 'Unknown')
                    }
                my_clinical: Clinical = Clinical(clinical_dict)
                the_clinical_info[my_clinical.bcr_patient_barcode] = my_clinical


def add_clinical_file(the_archive_file: str, the_internal_file: str, the_convert_type: str,
                      the_clinical_info: Dict[str, Clinical]) -> None:
    """
    Convert public clinical files into a dictionary of Clinical objects
    :param the_archive_file: pull path to directory with downloaded Public Clinical files
    :param the_internal_file: the name of the file inside the ZIP file to actually read
    :param the_convert_type: string describing file type to convert
    :param the_clinical_info: Dictionary of Clinical objects (key may be UUID or Barcode, depending on creation)
    :return:
    """
    print(f'add_clinical_file the_archive_file={the_archive_file} the_internal_file={the_internal_file} the_convert_type={the_convert_type}', flush=True)
    if 'convert-clinical_xml' == the_convert_type:
        add_clinical_file_clinical_xml(the_archive_file, the_internal_file, the_clinical_info)
    elif 'convert-clinical_json' == the_convert_type:
        add_clinical_file_clinical_json(the_archive_file, the_internal_file, the_clinical_info)
    elif 'convert-clinical_tsv' == the_convert_type:
        add_clinical_file_clinical_tsv(the_archive_file, the_internal_file, the_clinical_info)
    elif 'convert-clinical_xlsx' == the_convert_type:
        add_clinical_file_clinical_xlsx(the_archive_file, the_internal_file, the_clinical_info)
    else:
        add_error(f"ERROR add_clinical_file found unexpected convert type {the_convert_type}")


def write_headers_clinical(the_out_file: io.TextIOWrapper) -> None:
    """
    Write the headers for clinical objects
    Writes newline at the end.
    :param the_out_file: stream to which to write headers
    :return: None
    """
    write_tsv_list(the_out_file, CLINICAL_HEADERS, True, False)


def write_clinical_file(the_clinical_file: str, the_clinical_info: Dict[str, Clinical]) -> None:
    """
    Write dictionary of Clinical objects to a clinical/annotations TSV
    :param the_clinical_file: full path and filename to which to wrute clinical information
    :param the_clinical_info: Dictionary of Clinical objects
    :return: Nothing
    """
    print(f"write_clinical_file the_clinical_file={the_clinical_file}", flush=True)
    print(f"clinical sample count = {len(the_clinical_info)} file {the_clinical_file}", flush=True)
    if not len(the_clinical_info) > 2:
        print(f"ERROR generating clinical too few samples {len(the_clinical_info)} file {the_clinical_file}", flush=True)
    else:
        if not os.path.exists(os.path.dirname(the_clinical_file)):
            os.makedirs(os.path.dirname(the_clinical_file), exist_ok=True)
        clinical_list: List[Clinical] = list(the_clinical_info.values())
        clinical_list.sort(key=lambda my_clinical: my_clinical.bcr_patient_barcode, reverse=False)
        out_file: io.TextIOWrapper
        with open(the_clinical_file, 'w', encoding='utf-8') as out_file:
            write_headers_clinical(out_file)
            my_entry: Clinical
            for my_entry in clinical_list:
                my_entry.write_clinical(out_file)


def convert_clinical_batches(the_download_dir: str, the_clinical_file: str, the_clinical_files: List[GdcApiClinical]) -> None:
    """
    Take a list of GdcApiClinical files, with data downloaded, and write to the the_clinical_file TSV file
    :param the_download_dir: full path to directory with downloads
    :param the_clinical_file: full path and filename of TSV to write
    :param the_clinical_files: Python List of GdcApiClinical (file object) to write to TSV
    :return: Nothing
    """
    clinical_info: Dict[str, Clinical] = {}
    index: int
    clinical_file: GdcApiClinical
    for index, clinical_file in enumerate(the_clinical_files):
        if 0 == index % 100:
            print(f"convert_clinical_batches index={index}", flush=True)
        add_clinical_file(clinical_file.get_file_archive(the_download_dir),
                          clinical_file.file_name, clinical_file.get_convertable_type(),
                          clinical_info)
    write_clinical_file(the_clinical_file, clinical_info)
