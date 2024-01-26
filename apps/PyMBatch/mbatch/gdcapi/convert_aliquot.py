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
from typing import Dict, List, Tuple
import xmltodict
import pandas
from mbatch.test.common import xml_text, write_tsv_list, get_entry_as_list
from mbatch.test.common import get_dataframe_entry, add_warnings, remove_unsafe_characters
from mbatch.gdcapi.download_biospecimen import GdcApiBiospecimen
from mbatch.gdcapi.download_datafile import GdcApiDatafile, read_datafile_for_uuid_barcodes


# Sample	Patient	Type	BatchId	BCR	TSS	PlateId	AliquotCenter	ShipDate	SourceCenter
BATCH_HEADERS: List[str] = [
    'aliquot_barcode', 'aliquot_uuid',
    'patient_barcode', 'patient_uuid',
    'bcr', 'batch_id',
    'tissue_source_site', 'sex',
    'sample_type_name', 'sample_type_id',
    'ffpe', 'plate_id',
    'center_id', 'ship_date',
    'source_center'
]


# pylint: disable=too-many-instance-attributes,too-few-public-methods
class Aliquot:
    """
    Class to represent information about a single aliquot
    Samples are now called aliquots. An aliquot/sample is the most detailed barcode
    description of a sample.
    """
    # declare but do not set member attributes
    # use admin entries
    bcr: str
    batch_id: str
    # use patient entries
    patient_barcode: str
    patient_uuid: str
    tissue_source_site: str
    sex: str
    # use patient -> samples -> sample
    sample_type_name: str
    sample_type_id: str
    ffpe: str
    # use patient -> samples -> sample -> portions -> portion -> analytes -> analyte -> aliquots -> aliquot
    plate_id: str
    center_id: str
    ship_date: str
    source_center: str
    # from above source, and are the unique identifiers for this aliquot/sample
    aliquot_barcode: str
    aliquot_uuid: str

    # pylint: disable=too-many-arguments,too-many-locals
    def __init__(self: 'Aliquot', the_bcr: str, the_batch_id: str, the_patient_barcode: str,
                 the_patient_uuid: str, the_tissue_source_site: str, the_sex: str,
                 the_sample_type_name: str, the_sample_type_id: str, the_ffpe: str,
                 the_plate_id: str, the_center_id: str, the_ship_date: str,
                 the_source_center: str, the_aliquot_barcode: str, the_aliquot_uuid: str) -> None:
        """
        build an aliquot object. Use 'Unknown' in place of empty string.
        :param the_bcr: rough source of value documented in attribute declaration of class
        :param the_batch_id: rough source of value documented in attribute declaration of class
        :param the_patient_barcode: rough source of value documented in attribute declaration of class
        :param the_patient_uuid: rough source of value documented in attribute declaration of class
        :param the_tissue_source_site: rough source of value documented in attribute declaration of class
        :param the_sex: rough source of value documented in attribute declaration of class
        :param the_sample_type_name: rough source of value documented in attribute declaration of class
        :param the_sample_type_id: rough source of value documented in attribute declaration of class
        :param the_ffpe: rough source of value documented in attribute declaration of class
        :param the_plate_id: rough source of value documented in attribute declaration of class
        :param the_center_id: rough source of value documented in attribute declaration of class
        :param the_ship_date: rough source of value documented in attribute declaration of class
        :param the_source_center: rough source of value documented in attribute declaration of class
        :param the_aliquot_barcode: rough source of value documented in attribute declaration of class
        :param the_aliquot_uuid: rough source of value documented in attribute declaration of class
        """
        super().__init__()
        the_bcr = remove_unsafe_characters(the_bcr)
        the_batch_id = remove_unsafe_characters(the_batch_id)
        the_patient_barcode = remove_unsafe_characters(the_patient_barcode)
        the_patient_uuid = remove_unsafe_characters(the_patient_uuid)
        the_tissue_source_site = remove_unsafe_characters(the_tissue_source_site)
        the_sex = remove_unsafe_characters(the_sex)
        the_sample_type_name = remove_unsafe_characters(the_sample_type_name)
        the_sample_type_id = remove_unsafe_characters(the_sample_type_id)
        the_ffpe = remove_unsafe_characters(the_ffpe)
        the_plate_id = remove_unsafe_characters(the_plate_id)
        the_center_id = remove_unsafe_characters(the_center_id)
        the_ship_date = remove_unsafe_characters(the_ship_date)
        the_source_center = remove_unsafe_characters(the_source_center)
        the_aliquot_barcode = remove_unsafe_characters(the_aliquot_barcode)
        the_aliquot_uuid = remove_unsafe_characters(the_aliquot_uuid)
        # use admin entries
        self.bcr: str = 'Unknown' if '' == the_bcr else the_bcr
        self.batch_id: str = 'Unknown' if '' == the_batch_id else the_batch_id
        # use patient entries
        self.patient_barcode: str = 'Unknown' if '' == the_patient_barcode else the_patient_barcode
        self.patient_uuid: str = 'Unknown' if '' == the_patient_uuid else the_patient_uuid
        self.tissue_source_site: str = 'Unknown' if '' == the_tissue_source_site else the_tissue_source_site
        self.sex: str = 'Unknown' if '' == the_sex else the_sex
        # use patient -> samples -> sample
        self.sample_type_name: str = 'Unknown' if '' == the_sample_type_name else the_sample_type_name
        self.sample_type_id: str = 'Unknown' if '' == the_sample_type_id else the_sample_type_id
        self.ffpe: str = 'Unknown' if '' == the_ffpe else the_ffpe
        # use patient -> samples -> sample -> portions -> portion -> analytes -> analyte -> aliquots -> aliquot
        self.plate_id: str = 'Unknown' if '' == the_plate_id else the_plate_id
        self.center_id: str = 'Unknown' if '' == the_center_id else the_center_id
        self.ship_date: str = 'Unknown' if '' == the_ship_date else the_ship_date
        self.source_center: str = 'Unknown' if '' == the_source_center else the_source_center
        # from above source, and are the unique identifiers for this aliquot/sample
        self.aliquot_barcode: str = 'Unknown' if '' == the_aliquot_barcode else the_aliquot_barcode
        self.aliquot_uuid: str = 'Unknown' if '' == the_aliquot_uuid else the_aliquot_uuid
    # pylint: enable=too-many-arguments,too-many-locals

    def write_aliquot(self: 'Aliquot', the_out_file: io.TextIOWrapper) -> None:
        """
        Write out just the clinical portion of a file element to a file,
        after first calling parent history to write out its portion.
        Writes a newline.
        :param the_out_file: stream to which to write strings
        :return: None
        """
        data_list: List[str] = []
        for attr in BATCH_HEADERS:
            data_list.append(getattr(self, attr))
        write_tsv_list(the_out_file, data_list, True, False)
# pylint: enable=too-many-instance-attributes,too-few-public-methods

# ########################################################
# writing files
# ########################################################


# pylint: disable=too-many-nested-blocks,too-many-branches,too-many-locals,too-many-statements
def add_batch_bcr_xml(the_archive_file: str, the_internal_file: str, the_batch_info: Dict[str, Aliquot]) -> None:
    """
    Build batch information from BCR XML file. TCGA uses this.
    :param the_archive_file: the full path and filename for archive (ZIP) file
    :param the_internal_file: the name of the file inside the ZIP file to actually read
    :param the_batch_info: Dict of Aliquots (key may be UUID or Barcode, depending on creation)
    :return: Nothing
    """
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    with zipfile.ZipFile(the_archive_file, 'r') as zip_file:
        with zip_file.open(the_internal_file, mode="r") as in_file:
            # need type
            doc: Dict = xmltodict.parse(in_file.read())
            if 'ssf:tcga_bcr' in doc:
                print("add_biospecimen_file skip ssf file", flush=True)
            else:
                # use admin entries
                bcr: str = xml_text(doc['bio:tcga_bcr']['admin:admin']['admin:bcr'])
                batch_id: str = xml_text(doc['bio:tcga_bcr']['admin:admin']['admin:batch_number'])
                # use patient entries
                patient: Dict = doc['bio:tcga_bcr']['bio:patient']
                patient_barcode: str = xml_text(patient['shared:bcr_patient_barcode'])
                patient_uuid: str = xml_text(patient['shared:bcr_patient_uuid'])
                tissue_source_site: str = xml_text(patient['shared:tissue_source_site'])
                sex: str = xml_text(patient['shared:gender'])
                # use patient -> samples -> sample
                # sometimes bio_model:samples and sometimes bio:samples
                samples: List[Dict]
                # bio_model are not real samples
                # if 'bio_model:sample' in patient['bio:samples']:
                #    samples = get_entry_as_list(patient, 'bio:samples', 'bio_model:sample')
                # else:
                samples = get_entry_as_list(patient, 'bio:samples', 'bio:sample')
                samples_tmp: List[Dict] = get_entry_as_list(patient, 'bio:samples', 'bio_model:sample')
                samples.extend(samples_tmp)
                if len(samples) < 1:
                    add_warnings(f"No samples found in {the_archive_file} / {the_internal_file}")
                my_sample: Dict
                for my_sample in samples:
                    if isinstance(my_sample, dict):
                        sample_type_name: str = xml_text(my_sample['bio:sample_type'])
                        sample_type_id: str = xml_text(my_sample['bio:sample_type_id'])
                        analytes: List[Dict] = []
                        if 'bio:portions' in my_sample:
                            # add bio:portions -> bio:shipment_portion to aliquot list directly (RPPA samples are marked as this)
                            ship_portions: List[Dict] = get_entry_as_list(my_sample, 'bio:portions', 'bio:shipment_portion')
                            my_ship_portion: Dict
                            for my_ship_portion in ship_portions:
                                plate_id: str = xml_text(my_ship_portion['bio:plate_id'])
                                center_id: str = xml_text(my_ship_portion['bio:center_id'])
                                day: str = xml_text(my_ship_portion['bio:shipment_portion_day_of_shipment'])
                                month: str = xml_text(my_ship_portion['bio:shipment_portion_month_of_shipment'])
                                year: str = xml_text(my_ship_portion['bio:shipment_portion_year_of_shipment'])
                                ship_date: str = 'Unknown'
                                if isinstance(year, str):
                                    ship_date = year + "-" + month + "-" + day
                                source_center: str = "Unknown"
                                aliquot_barcode: str = xml_text(my_ship_portion['bio:shipment_portion_bcr_aliquot_barcode'])
                                aliquot_uuid: str = xml_text(my_ship_portion['bio:bcr_shipment_portion_uuid'])
                                ffpe: str = xml_text(my_ship_portion['bio:is_ffpe'])
                                b_batch: Aliquot = Aliquot(bcr, batch_id, patient_barcode,
                                                           patient_uuid, tissue_source_site, sex,
                                                           sample_type_name, sample_type_id, ffpe,
                                                           plate_id, center_id, ship_date,
                                                           source_center, aliquot_barcode, aliquot_uuid)
                                the_batch_info[b_batch.aliquot_uuid] = b_batch
                            # add bio:portions -> bio:portion to analyte list (other samples are marked as this)
                            portions: List[Dict] = get_entry_as_list(my_sample, 'bio:portions', 'bio:portion')
                            my_portion: Dict
                            for my_portion in portions:
                                tmp_analytes: List[Dict] = get_entry_as_list(my_portion, 'bio:analytes', 'bio:analyte')
                                tmp_analyte: Dict
                                for tmp_analyte in tmp_analytes:
                                    analytes.append(tmp_analyte)
                        else:
                            analytes = get_entry_as_list(my_sample, 'bio:analytes', 'bio:analyte')
                        my_analyte: Dict
                        for my_analyte in analytes:
                            # use analyte -> aliquots -> aliquot
                            if isinstance(my_analyte, dict):
                                aliquots: List[Dict] = get_entry_as_list(my_analyte, 'bio:aliquots', 'bio:aliquot')
                                my_aliquot: Dict
                                for my_aliquot in aliquots:
                                    if isinstance(my_aliquot, dict):
                                        plate_id: str = xml_text(my_aliquot['bio:plate_id'])
                                        center_id: str = xml_text(my_aliquot['bio:center_id'])
                                        day: str = xml_text(my_aliquot['bio:day_of_shipment'])
                                        month: str = xml_text(my_aliquot['bio:month_of_shipment'])
                                        year: str = xml_text(my_aliquot['bio:year_of_shipment'])
                                        ship_date: str = 'Unknown'
                                        if isinstance(year, str):
                                            ship_date = year + "-" + month + "-" + day
                                        source_center: str = xml_text(my_aliquot['bio:source_center'])
                                        aliquot_barcode: str = xml_text(my_aliquot['bio:bcr_aliquot_barcode'])
                                        aliquot_uuid: str = xml_text(my_aliquot['bio:bcr_aliquot_uuid'])
                                        ffpe: str = xml_text(my_aliquot['bio:is_derived_from_ffpe'])
                                        b_batch: Aliquot = Aliquot(bcr, batch_id, patient_barcode,
                                                                   patient_uuid, tissue_source_site, sex,
                                                                   sample_type_name, sample_type_id, ffpe,
                                                                   plate_id, center_id, ship_date,
                                                                   source_center, aliquot_barcode, aliquot_uuid)
                                        the_batch_info[b_batch.aliquot_uuid] = b_batch
# pylint: enable=too-many-nested-blocks,too-many-branches,too-many-locals,too-many-statements


# pylint: disable=too-many-nested-blocks,too-many-branches,too-many-locals,too-many-statements
def add_batch_cgci_xlsx(the_archive_file: str, the_internal_file: str, the_batch_info: Dict[str, Aliquot],
                         the_datafile_list: List[GdcApiDatafile],
                         the_sample_dir: str) -> None:
    """
    Build batch information from XLSX files. CGCI Program uses this.
    CGCI only puts barcodes in their XLSX files, so we have to pass in map of barcodes to UUIDs.
    :param the_archive_file: the full path and filename for archive (ZIP) file
    :param the_internal_file: the name of the file inside the ZIP file to actually read
    :param the_batch_info: Dict of Aliquots (key may be UUID or Barcode, depending on creation)
    :param the_datafile_list: List of GdcApiDatafile objects (REQUIRED for some conversions)
    :param the_sample_dir: full path to sample dir, which are associated with the dataset index (REQUIRED for some conversions)
    :return: Nothing
    """
    # read barcode to uuid conversion data
    barcode_to_uuid: Dict[str, str]
    (_, barcode_to_uuid) = read_datafile_for_uuid_barcodes(the_datafile_list, the_sample_dir)
    # convert given file
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    with zipfile.ZipFile(the_archive_file, 'r') as zip_file:
        with zip_file.open(the_internal_file, mode="r") as in_file:
            # add ,0 to read_excel to resolve PyLint mapping issue
            my_df: pandas.DataFrame = pandas.read_excel(in_file, 0, dtype='str')
            # convert missing values to empty string
            my_df = my_df.mask(my_df.isna(), '')
            # iterate rows of dataframe
            index: int
            row: pandas.Series
            for index, row in my_df.iterrows():
                if 0 == index % 1000:
                    print(f"add_batch_cgci_xlsx index={index}", flush=True)
                bcr: str = get_dataframe_entry(row, '', 'Unknown')
                batch_id: str = get_dataframe_entry(row, '', 'Unknown')
                patient_barcode: str = get_dataframe_entry(row, 'Case ID', 'Unknown')
                # get UUID from dataset-samples entries
                patient_uuid: str = barcode_to_uuid.get(patient_barcode, 'Unknown')
                tissue_source_site: str = get_dataframe_entry(row, '', 'Unknown')
                sex: str = get_dataframe_entry(row, '', 'Unknown')
                percent_tumor: str = get_dataframe_entry(row, '% tumour top', 'Unknown')
                # do not initialize, is set in if/else below
                sample_type_name: str
                sample_type_id: str
                if '0' == percent_tumor:
                    sample_type_name = 'normal'
                    sample_type_id = 'normal'
                else:
                    sample_type_name = 'tumor'
                    sample_type_id = 'tumor'
                ffpe: str = get_dataframe_entry(row, 'is_ffpe', 'Unknown')
                plate_id: str = get_dataframe_entry(row, 'Shipment ID', 'Unknown')
                center_id: str = get_dataframe_entry(row, 'Extraction method', 'Unknown')
                ship_date: str = get_dataframe_entry(row, '', 'Unknown')
                source_center: str = get_dataframe_entry(row, '', 'Unknown')
                aliquot_barcode: str = get_dataframe_entry(row, 'Aliquot ID', 'Unknown')
                # get UUID from dataset-samples entries
                aliquot_uuid: str = barcode_to_uuid.get(aliquot_barcode, 'Unknown')
                b_batch: Aliquot = Aliquot(bcr, batch_id, patient_barcode,
                                           patient_uuid, tissue_source_site, sex,
                                           sample_type_name, sample_type_id, ffpe,
                                           plate_id, center_id, ship_date,
                                           source_center, aliquot_barcode, aliquot_uuid)
                the_batch_info[b_batch.aliquot_uuid] = b_batch
# pylint: enable=too-many-nested-blocks,too-many-branches,too-many-locals,too-many-statements


# pylint: disable=too-many-nested-blocks,too-many-branches,too-many-locals,too-many-statements
def add_batch_fm_tsv(the_archive_file: str, the_internal_file: str, the_batch_info: Dict[str, Aliquot]) -> None:
    """
    Build batch information from TSV files. FM Program uses this..
    :param the_archive_file: the full path and filename for archive (ZIP) file
    :param the_internal_file: the name of the file inside the ZIP file to actually read
    :param the_batch_info: Dict of Aliquots (key may be UUID or Barcode, depending on creation)
    :return: Nothing
    """
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    with zipfile.ZipFile(the_archive_file, 'r') as zip_file:
        with zip_file.open(the_internal_file, mode="r") as in_file:
            index: int = 0
            keys: List[str]
            for line in io.TextIOWrapper(in_file, encoding="utf-8"):
                if 0 == index % 1000:
                    print(f"add_batch_fm_tsv index={index}", flush=True)
                if 0 == index:
                    keys = line.rstrip('\n').split('\t')
                else:
                    values: List[str] = line.rstrip('\n').split("\t")
                    batch_line: Dict[str, str] = dict(zip(keys, values))
                    bcr: str = batch_line.get('', 'Unknown')
                    batch_id: str = batch_line.get('', 'Unknown')
                    patient_barcode: str = batch_line.get('cases.submitter_id', 'Unknown')
                    patient_uuid: str = batch_line.get('case_id', 'Unknown')
                    tissue_source_site: str = batch_line.get('cases.primary_site', 'Unknown')
                    sex: str = batch_line.get('', 'Unknown')
                    sample_type_name: str = batch_line.get('samples.tumor_descriptor', 'Unknown')
                    sample_type_id: str = batch_line.get('samples.sample_type', 'Unknown')
                    ffpe: str = batch_line.get('', 'Unknown')
                    plate_id: str = batch_line.get('read_groups.library_strategy', 'Unknown')
                    center_id: str = batch_line.get('read_groups.platform', 'Unknown')
                    ship_date: str = batch_line.get('', 'Unknown')
                    source_center: str = batch_line.get('read_groups.sequencing_center', 'Unknown')
                    aliquot_barcode: str = batch_line.get('aliquots.submitter_id', 'Unknown')
                    aliquot_uuid: str = batch_line.get('aliquot_id', 'Unknown')
                    b_batch: Aliquot = Aliquot(bcr, batch_id, patient_barcode,
                                               patient_uuid, tissue_source_site, sex,
                                               sample_type_name, sample_type_id, ffpe,
                                               plate_id, center_id, ship_date,
                                               source_center, aliquot_barcode, aliquot_uuid)
                    the_batch_info[b_batch.aliquot_uuid] = b_batch
                index += 1
# pylint: enable=too-many-nested-blocks,too-many-branches,too-many-locals,too-many-statements


def build_target_barcode_tuples(the_row: pandas.Series, the_key: str, the_type: str, the_barcode_type: List[Tuple[str, str]]) -> None:
    """
    Build barcode to type tuples for use by TARGET conversion
    :param the_row: Pandas Series (row from file) being converted to an aliquot
    :param the_key: Key (column name) to get from row
    :param the_type: Type (from elsewhere) for this Barcode.
    :param the_barcode_type: List of barcodes and their types
    :return: Nothing
    """
    entry: str = get_dataframe_entry(the_row, the_key, '')
    if '' != entry:
        barcode: str
        for barcode in entry.split(','):
            the_barcode_type.append((barcode, the_type))


# pylint: disable=too-many-nested-blocks,too-many-branches,too-many-locals,too-many-statements
def add_batch_target_xlsx(the_archive_file: str, the_internal_file: str, the_batch_info: Dict[str, Aliquot],
                          the_datafile_list: List[GdcApiDatafile], the_sample_dir: str) -> None:
    """
    Build batch information from XLSX files. TARGET Program uses this.
    TARGET only puts barcodes in their XLSX files, so we have to pass in map of barcodes to UUIDs.
    :param the_archive_file: the full path and filename for archive (ZIP) file
    :param the_internal_file: the name of the file inside the ZIP file to actually read
    :param the_batch_info: Dict of Aliquots (key may be UUID or Barcode, depending on creation)
    :param the_datafile_list: List of GdcApiDatafile objects (REQUIRED for some conversions)
    :param the_sample_dir: full path to sample dir, which are associated with the dataset index (REQUIRED for some conversions)
    :return: Nothing
    """
    # read barcode to uuid conversion data
    barcode_to_uuid: Dict[str, str]
    (_, barcode_to_uuid) = read_datafile_for_uuid_barcodes(the_datafile_list, the_sample_dir)
    # convert given file
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    with zipfile.ZipFile(the_archive_file, 'r') as zip_file:
        with zip_file.open(the_internal_file, mode="r") as in_file:
            # add ,0 to read_excel to resolve PyLint mapping issue
            my_df: pandas.DataFrame = pandas.read_excel(in_file, 0)
            # convert missing values to empty string
            my_df = my_df.mask(my_df.isna(), '')
            # iterate rows of dataframe
            index: int
            row: pandas.Series
            for index, row in my_df.iterrows():
                if 0 == index % 1000:
                    print(f"add_batch_target_xlsx index={index}", flush=True)
                # skip first row -- due to double row of headers
                if index > 1:
                    bcr: str = get_dataframe_entry(row, '', 'Unknown')
                    batch_id: str = get_dataframe_entry(row, '', 'Unknown')
                    patient_barcode: str = get_dataframe_entry(row, 'Case USI', 'Unknown')
                    # get UUID from dataset-samples entries
                    patient_uuid: str = barcode_to_uuid.get(patient_barcode, 'Unknown')
                    tissue_source_site: str = get_dataframe_entry(row, '', 'Unknown')
                    sex: str = get_dataframe_entry(row, '', 'Unknown')
                    sample_type_name: str = get_dataframe_entry(row, 'Comments', 'Unknown')
                    ffpe: str = get_dataframe_entry(row, '', 'Unknown')
                    # incorrect spelling from spreadsheets
                    plate_id: str = get_dataframe_entry(row, 'Multiple RNA Aliqots From the Same Tissue (Y/N)', 'Unknown')
                    center_id: str = get_dataframe_entry(row, 'Multiple DNA Aliquots From the Same Tissue (Y/N)', 'Unknown')
                    ship_date: str = get_dataframe_entry(row, '', 'Unknown')
                    source_center: str = get_dataframe_entry(row, '', 'Unknown')
                    # collect barcode type tuples
                    barcode_type: List[(str, str)] = []
                    build_target_barcode_tuples(row, 'Diagnostic Tumor RNA Sample ID', 'RNA Tumor', barcode_type)
                    build_target_barcode_tuples(row, 'Matched Normal RNA Sample ID', 'RNA Normal', barcode_type)
                    build_target_barcode_tuples(row, 'Relapse Tumor RNA Sample ID', 'RNA Tumor Relapse', barcode_type)
                    build_target_barcode_tuples(row, 'Diagnostic Tumor DNA Sample ID', 'DNA Tumor', barcode_type)
                    build_target_barcode_tuples(row, 'Matched Normal DNA Sample ID', 'DNA Normal', barcode_type)
                    build_target_barcode_tuples(row, 'Relapse Tumor DNA Sample ID', 'DNA Tumor Relapse', barcode_type)
                    aliquot_barcode: str
                    sample_type_id: str
                    for aliquot_barcode, sample_type_id in barcode_type:
                        # get UUID from dataset-samples entries
                        aliquot_uuid: str = barcode_to_uuid.get(aliquot_barcode, 'Unknown')
                        b_batch: Aliquot = Aliquot(bcr, batch_id, patient_barcode,
                                                   patient_uuid, tissue_source_site, sex,
                                                   sample_type_name, sample_type_id, ffpe,
                                                   plate_id, center_id, ship_date,
                                                   source_center, aliquot_barcode, aliquot_uuid)
                        the_batch_info[b_batch.aliquot_barcode] = b_batch
# pylint: enable=too-many-nested-blocks,too-many-branches,too-many-locals,too-many-statements


# pylint: disable=too-many-arguments
def add_biospecimen_file(the_archive_file: str, the_internal_file: str, the_convert_type: str,
                         the_batch_info: Dict[str, Aliquot],
                         the_datafile_list: List[GdcApiDatafile],
                         the_sample_dir: str) -> None:
    """
    Add aliquot(s) from passed in biospecimen file into passed in the_batch_info Dictionary
    :param the_archive_file: the full path and filename for archive (ZIP) file
    :param the_internal_file: the name of the file inside the ZIP file to actually read
    :param the_convert_type: string giving type of convert to do
    :param the_batch_info: Dict of Aliquots (key may be UUID or Barcode, depending on creation)
    :param the_datafile_list: List of GdcApiDatafile objects (REQUIRED for some conversions)
    :param the_sample_dir: full path to sample dir, which are associated with the dataset index (REQUIRED for some conversions)
    :return: Nothing
    """
    print(f'add_biospecimen_file the_archive_file={the_archive_file} the_internal_file={the_internal_file} the_convert_type={the_convert_type}', flush=True)
    if 'convert-bcr_xml' == the_convert_type:
        add_batch_bcr_xml(the_archive_file, the_internal_file, the_batch_info)
    elif 'convert-cgci_xlsx' == the_convert_type:
        add_batch_cgci_xlsx(the_archive_file, the_internal_file, the_batch_info, the_datafile_list, the_sample_dir)
    elif 'convert-fm_tsv' == the_convert_type:
        add_batch_fm_tsv(the_archive_file, the_internal_file, the_batch_info)
    elif 'convert-target_xlsx' == the_convert_type:
        add_batch_target_xlsx(the_archive_file, the_internal_file, the_batch_info, the_datafile_list, the_sample_dir)
# pylint: enable=too-many-arguments


def write_headers_batch(the_out_file: io.TextIOWrapper) -> None:
    """
    Write the headers for batch objects
    Writes newline at the end.
    :param the_out_file: stream to which to write headers
    :return: None
    """
    write_tsv_list(the_out_file, BATCH_HEADERS, True, False)


def write_batches_file(the_batch_file: str, the_batch_info: Dict[str, Aliquot]) -> None:
    """
    Write a batches.tsv file out from the biospecimen built Dictionary of Aliquots.
    This is a program/project level batches.tsv that will get trimmed for a dataset.
    :param the_batch_file: full path and filename to which to write
    :param the_batch_info: Dictionary of Aliquots to write a batch file
    :return: Nothing
    """
    print(f"write_batches_file the_batch_file={the_batch_file}", flush=True)
    if not len(the_batch_info) > 10:
        print(f"ERROR generating batches too few samples {len(the_batch_info)} file {the_batch_file}", flush=True)
    else:
        if not os.path.exists(os.path.dirname(the_batch_file)):
            os.makedirs(os.path.dirname(the_batch_file), exist_ok=True)
        out_file: io.TextIOWrapper
        aliquot_list: List[Aliquot] = list(the_batch_info.values())
        aliquot_list.sort(key=lambda my_aliquot: my_aliquot.aliquot_barcode, reverse=False)
        with open(the_batch_file, 'w', encoding='utf-8') as out_file:
            write_headers_batch(out_file)
            my_entry: Aliquot
            for my_entry in aliquot_list:
                my_entry.write_aliquot(out_file)


def convert_biospecimen_batches(the_download_dir: str, the_batch_file: str,
                                the_biospecimen_files: List[GdcApiBiospecimen],
                                the_datafile_list: List[GdcApiDatafile],
                                the_sample_dir: str) -> None:
    """
    Convert biospecimen files into a dictionary of GdcApiBiospecimen objects
    :param the_download_dir: pull path to directory with downloaded biospecimen files
    :param the_batch_file: full path and filename of program/project level batches.tsv to write.
    :param the_biospecimen_files: List of GdcApiBiospecimen objects
    :param the_datafile_list: List of GdcApiDatafile objects (REQUIRED for some conversions)
    :param the_sample_dir: full path to sample dir, which are associated with the dataset index (REQUIRED for some conversions)
    :return:
    """
    batch_info: Dict[str, Aliquot] = {}
    index: int
    biospecimen_file: GdcApiBiospecimen
    for index, biospecimen_file in enumerate(the_biospecimen_files):
        if 0 == index % 100:
            print(f"convert_biospecimen_batches index={index}", flush=True)
        add_biospecimen_file(biospecimen_file.get_file_archive(the_download_dir),
                             biospecimen_file.file_name, biospecimen_file.get_convertable_type(),
                             batch_info, the_datafile_list, the_sample_dir)
    write_batches_file(the_batch_file, batch_info)
