#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright (c) 2011-2022 University of Texas MD Anderson Cancer Center

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


from pathlib import Path
from typing import List, Dict
import typing
import math
import hashlib
import io
import os
import re
import shutil
import datetime
import zipfile
import pandas


def timestamp_file_if_exists(the_file_path: str, the_timestamp: str) -> None:
    """
    if file exists, rename with the given timestamp
    :param the_file_path: Check if this file exists
    :param the_timestamp: if it exists, rename with _<the_timestamp> on the end
    :return: nothing
    """
    if os.path.exists(the_file_path):
        os.rename(the_file_path, f"{the_file_path}_{the_timestamp}")


def generate_file_md5(the_file_path: str, the_block_size=2**20) -> str:
    """
    based on https://stackoverflow.com/questions/1131220/get-md5-hash-of-big-files-in-python
    :param the_file_path: full path to file to get md5
    :param the_block_size: block size to read at one time
    :return: string MD5 for file
    """
    md5lib: hashlib = hashlib.md5()
    my_file: io.BufferedReader
    with open(os.path.join(the_file_path), "rb") as my_file:
        while True:
            buf: bytes = my_file.read(the_block_size)
            if not buf:
                break
            md5lib.update(buf)
    return md5lib.hexdigest()


def read_file_to_string(the_file: str) -> str:
    """
    Open and read file into a string -- strips whitespace from front and back
    :param the_file: file to read
    :return: string from file -- empty if file does not exist
    """
    file_contents: str
    if Path(the_file).exists():
        with open(the_file, 'r', encoding='utf-8') as my_file:
            file_contents = my_file.read()
    else:
        file_contents = ""
    file_contents = file_contents.strip()
    return file_contents


# pylint: disable=too-many-branches
def get_filenames(the_dir: str) -> str:
    """
    Get sorted, unique list of filenames in directory structure
    :param the_dir: directory to check for filenames
    :return: list of filenames, sorted and unique
    """
    matching_files: List[str] = []
    for _, _, file_names in os.walk(the_dir):
        check_name: str
        for check_name in file_names:
            if check_name.endswith("ngchm"):
                matching_files.append("NGCHM")
            elif check_name == "HCData.tsv":
                matching_files.append("HierarchicalClustering")
            elif check_name.startswith("BoxPlot"):
                matching_files.append("BoxPlot")
            elif check_name == "CDP_Plot_Data1_Diagram.PNG":
                matching_files.append("CDP")
            elif check_name == "DSCOverview.tsv":
                matching_files.append("DSCOverview.tsv")
            elif check_name.startswith("SupervisedClust"):
                matching_files.append("SupervisedClustering")
            elif check_name == "KW_Dunns_Diagram.tsv":
                matching_files.append("Kruskal-Wallis Dunn's Test")
            elif check_name == "PCAValues.tsv":
                matching_files.append("PCA")
            elif check_name.startswith("UMAP_Data"):
                matching_files.append("UMAP")
    if len(matching_files) > 0:
        matching_files = list(set(matching_files))
        matching_files.sort()
    result: str = ""
    for check_name in matching_files:
        if "" == result:
            result = result + check_name
        else:
            result = result + "|" + check_name
    return result
# pylint: enable=too-many-branches


def get_newest_file_list(the_dir: str, the_file: str) -> List[str]:
    """
    Find the files (sorted in order) in given directory structure.
    Pick ones with same last directory in path.
    Directories will be named such that reverse sort, Z or 9 first, gives newest.
    Return full path.
    :param the_dir: directory to check for newest file
    :param the_file: filename to search for
    :return: full path to newest file
    """
    newest_files: List[str] = []
    matching_files: List[str] = []
    current_dir: str
    file_names: List[str]
    for current_dir, _, file_names in os.walk(the_dir):
        if the_file in file_names:
            matching_files.append(f"{current_dir}/{the_file}")
    if len(matching_files) > 0:
        matching_files.sort(reverse=True)
        current_file: str
        comp_file: str = ""
        for current_file in matching_files:
            if "" == comp_file:
                index_a: int = current_file.rfind("/")
                tmp: str = current_file[:index_a]
                index_b: int = tmp.rfind("/")
                tmp = current_file[:index_b]
                index_c: int = tmp.rfind("/")
                comp_file = current_file[index_c:]
                newest_files.append(current_file)
            else:
                if current_file.endswith(comp_file):
                    newest_files.append(current_file)
    return newest_files


def get_newest_file(the_dir: str, the_file: str) -> typing.Optional[str]:
    """
    Find the newest file (based on sort order) in given directory structure.
    Directories will be named such that reverse sort, Z or 9 first, gives newest.
    Return full path.
    :param the_dir: directory to check for newest file
    :param the_file: filename to search for
    :return: full path to newest file
    """
    newest_file: typing.Optional[str] = None
    matching_files: List[str] = []
    current_dir: str
    file_names: List[str]
    for current_dir, _, file_names in os.walk(the_dir):
        if the_file in file_names:
            matching_files.append(f"{current_dir}/{the_file}")
    if len(matching_files) > 0:
        matching_files.sort(reverse=True)
        newest_file = matching_files[0]
    return newest_file


def get_newest_dir(the_dir: str) -> typing.Optional[str]:
    """
    Find the newest directory (based on sort order) in given directory.
    Directories will be named such that reverse sort, Z or 9 first, gives newest.
    Return full path.
    :param the_dir: directory to check for newest subdirectory
    :return: full path to newest subdirectory
    """
    newest_dir: typing.Optional[str] = None
    if os.path.exists(the_dir):
        dir_entry: os.DirEntry
        subdir_list: List[str] = [dir_entry.path for dir_entry in os.scandir(the_dir) if dir_entry.is_dir()]
        subdir_list.sort(reverse=True)
        newest_dir = subdir_list[0]
    return newest_dir


def get_path(the_dir: os.DirEntry) -> str:
    """
    get path entry from class
    used to pass key function to array processing
    :param the_dir: os.DirEntry instance
    :return: path attribute value
    """
    return the_dir.path


def get_sorted_dirs(the_dir: str) -> List[os.DirEntry]:
    """
    Return sorted directory listing
    Return full path.
    :param the_dir: directory to check for newest subdirectory
    :return: sorted list of os.DirEntry elements
    """
    dir_entry: os.DirEntry
    subdir_list: List[os.DirEntry] = [dir_entry for dir_entry in os.scandir(the_dir) if dir_entry.is_dir()]
    subdir_list.sort(key=get_path, reverse=False)
    return subdir_list


def next_sub_dir_starts_with(the_dir: str, the_prefix: str) -> bool:
    """
    Determine if the first subdirectory (sorted) starts with the given prefix
    :param the_dir: string, full path to directory to check for subdirectories
    :param the_prefix: prefix to check for in subdirectories
    :return: True if first subdirectory starts with prefix. False otherwise
    """
    # print(f"next_sub_dir_starts_with the_dir={the_dir}", flush=True)
    # print(f"next_sub_dir_starts_with the_prefix={the_prefix}", flush=True)
    my_dirs: List[os.DirEntry] = get_sorted_dirs(the_dir)
    matches: bool = False
    if len(my_dirs) > 0:
        if my_dirs[0].name.startswith(the_prefix):
            matches = True
    return matches


def get_sorted_files(the_dir: str) -> List[os.DirEntry]:
    """
    Return sorted directory listing of files only
    Return full path.
    :param the_dir: directory to get list of files
    :return: sorted list of os.DirEntry elements which are files
    """
    dir_entry: os.DirEntry
    file_list: List[os.DirEntry] = [dir_entry for dir_entry in os.scandir(the_dir) if dir_entry.is_file()]
    file_list.sort(key=get_path, reverse=False)
    return file_list


def copy_dirs_and_files(the_source_dir: str, the_dest_dir: str) -> None:
    """
    Copy contents of source dir to the destination dir
    :param the_source_dir: full path of directory with contents to copy
    :param the_dest_dir: full path to directory to copy contents
    :return: None
    """
    # https://stackoverflow.com/questions/33790444/python-recursively-copy-batch-of-directories
    # print(f"copy_dirs_and_files the_source_dir={the_source_dir}", flush=True)
    # print(f"copy_dirs_and_files the_dest_dir={the_dest_dir}", flush=True)
    if not os.path.exists(the_dest_dir):
        os.makedirs(the_dest_dir)
    for candidate in os.listdir(the_source_dir):
        src_can = os.path.join(the_source_dir, candidate)
        dst_can = os.path.join(the_dest_dir, candidate)
        if os.path.isdir(src_can):
            copy_dirs_and_files(src_can, dst_can)
        else:
            shutil.copyfile(src_can, dst_can)


def get_current_timestamp() -> str:
    """
    YYYY-MM-DD-HHMM time format string for now
    :return: YYYY-MM-DD-HHMM time format string for now
    """
    return datetime.datetime.today().strftime('%Y-%m-%d-%H%M')


def write_tsv_list(the_out_file: io.TextIOWrapper, the_list: List[str], the_new_line: bool, the_trailing_tab: bool) -> None:
    """
    Write to a TSV file - can be a full or partial line
    :param the_out_file: stream to write to
    :param the_list: list of tokens/strings to write to line
    :param the_new_line: if True, write a newline on the end
    :param the_trailing_tab: if True, write a tab on the end
    :return: Nothing
    """
    index: int
    element: str
    for index, element in enumerate(the_list):
        if index > 0:
            the_out_file.write("\t")
        # everything written must be a string
        the_out_file.write(f'{element}')
    if the_new_line:
        the_out_file.write("\n")
    if the_trailing_tab:
        the_out_file.write("\t")


def read_headers(the_in_file: io.TextIOWrapper) -> List[str]:
    """
    read headers, splitting the line
    Be sure to string \n off end
    :param the_in_file: stream to read
    :return: List of strings/tokens from file
    """
    header_line: str = the_in_file.readline()
    if isinstance(header_line, bytes):
        header_line = header_line.decode('utf-8')
    header_line = header_line.rstrip('\n')
    return header_line.split("\t")


def xml_text(the_xml_element: Dict, the_default: typing.Optional[str] = '') -> str:
    """
    check if the element has a #text attribute, return value, or ''
    :param the_xml_element: a Dict created by xmltodict
    :param the_default: default string to use ('' if none is provided)
    :return: value of #text or the default value
    """
    return the_xml_element.get('#text', the_default)


def get_data_string_text(the_xml_element: Dict, the_attr: str, the_default) -> str:
    """
    This is used in combination with xml_text. Check to see if the_attr exists in the Dictionary object.
    If it does, get the value and pass the value and the default to xml_text. Otherwise, return default.
    xml_text checks for a #text element, or uses default.
    :param the_xml_element: a Dict created by xmltodict
    :param the_attr: attribute to get from the Dictionary object
    :param the_default: default value to use
    :return: value of #text of requested attribute, or the default value
    """
    val: str = the_default
    if the_attr in the_xml_element:
        val = xml_text(the_xml_element[the_attr], the_default)
    return val


def get_entry_as_list(the_parent: Dict, the_list_tag: str, the_element_tag: str) -> List[Dict]:
    """
    xmltodict and Dict have problems with lists with a single element. This function grabs the elements
    within the_list_tag named with the_element_tag, and returns them as a List of dictionary
    objects, even if there is only one element in the List.
    :param the_parent: dictionary object
    :param the_list_tag: tag wrapping the list elements
    :param the_element_tag: tag marking list elements
    :return: List of dictionary objects with given element tag
    """
    return_val: List[Dict] = []
    # not, no type given as type could be Dict or List[Dict]
    list_value = the_parent[the_list_tag]
    if list_value is not None:
        if the_element_tag in list_value:
            raw_value = list_value[the_element_tag]
            if isinstance(raw_value, list):
                return_val = raw_value
            else:
                return_val.append(raw_value)
    return return_val


def get_tsv_entry_dict(the_dict: Dict[str, str], the_key: str, the_default: str) -> str:
    """
    get the given attribute from the dictionary, using the default value if needed.
    Then if the attribute value was an empty string, use the default value instead.
    Also uses default if key is empty string.
    :param the_dict: dictionary with potentially a particular key
    :param the_key: attribute to use from dictionary
    :param the_default: default to use if no key or value is empty string
    :return: value (but not empty string) or default
    """
    value: str = the_default
    if '' != the_key:
        if the_key in the_dict:
            value = the_dict[the_key]
    if '' == value:
        value = the_default
    return value


def get_dataframe_entry(the_series: pandas.Series, the_key: str, the_default: str) -> str:
    """
    get the given attribute from the pandas.Series, using the default value if needed.
    Then if the attribute value was an empty string, use the default value instead.
    Also uses default if key is empty string.
    :param the_series: pandas.Series with potentially a particular key
    :param the_key: attribute to use from pandas.Series
    :param the_default: default to use if no key or value is empty string
    :return: value (but not empty string) or default
    """
    value: str = the_default
    if '' != the_key:
        if the_key in the_series:
            value = the_series[the_key]
    if '' == value:
        value = the_default
    return value


def handle_error(the_error) -> None:
    """
    print error
    :param the_error: passed from multiprocessing pool
    :return: None
    """
    print("=====================", flush=True)
    print(f"type(the_error)={type(the_error)}", flush=True)
    print(f"the_error={the_error}", flush=True)
    print("---------------------", flush=True)


WARN_LIST: List[str] = []
ERROR_LIST: List[str] = []


def print_warnings() -> None:
    """
    Print strings on WARN_LIST.
    :return: Nothing
    """
    my_str: str
    for my_str in WARN_LIST:
        print(my_str, flush=True)


def print_errors() -> None:
    """
    Print strings on ERROR_LIST.
    :return: Nothing
    """
    my_str: str
    for my_str in ERROR_LIST:
        print(my_str, flush=True)


def add_warnings(the_string: str) -> None:
    """
    Add given string to list of warnings.
    :param the_string: Add given string to list of warnings
    :return: nothing
    """
    global WARN_LIST
    print(the_string, flush=True)
    WARN_LIST.append(the_string)


def add_error(the_string: str) -> None:
    """
    Add given string to list of errors.
    :param the_string: Add given string to list of errors
    :return: nothing
    """
    global ERROR_LIST
    print(the_string, flush=True)
    ERROR_LIST.append(the_string)


def delete_from_dirs(the_root_dir: str, the_delete_filename: str) -> None:
    """
    Search root dir and delete all instances of given filename.
    Used mostly to clean up Rplot.pdf files left over from other packages.
    :param the_root_dir:
    :param the_delete_filename:
    :return: nothing
    """
    my_dir: str
    # _: List[str] ignored variable
    file_names: List[str]
    file_name: str
    for my_dir, _, file_names in os.walk(the_root_dir):
        for file_name in file_names:
            if file_name == the_delete_filename:
                print(f'Delete file {my_dir}{os.sep}{file_name}', flush=True)
                os.remove(f'{my_dir}{os.sep}{file_name}')


def delete_directory_contents(the_dir: str) -> None:
    """"
    Delete contents of directory.
    No simple way in Python, so I wrote the.
    :param the_dir: Delete contents of this directory, but not this directory.
    :return: Nothing
    """
    filename: str
    if os.path.exists(the_dir):
        for filename in os.listdir(the_dir):
            file_path: str = os.path.join(the_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)


def extract_zip_to_dir(the_zip_path: str, the_out_dir: str) -> None:
    """
    Extract the given ZIP to the given directory
    :param the_zip_path: ZIP file to extract
    :param the_out_dir: Directory to which to extract
    :return: Nothing
    """
    zip_file: zipfile.ZipFile
    with zipfile.ZipFile(the_zip_path, 'r') as zip_file:
        zip_file.extractall(the_out_dir)


def archive_dir_contents_to_zip(the_zip_path: str, the_archive_dir: str) -> None:
    """
    Compress contents of directory into ZIP
    :param the_zip_path: ZIP file to create
    :param the_archive_dir: Compress contents of this directory
    :return: nothing
    """
    # file comes in with a .zip on the end
    no_zip_path: str = the_zip_path[:-4]
    shutil.make_archive(no_zip_path, 'zip', the_archive_dir)


def copy_archive_file_to_regular(the_zip_path: str, the_internal_file: str, the_external_file: str) -> None:
    """
    Copy a file in the ZIP file to the outside
    :param the_zip_path: Path to ZIP file (with name)
    :param the_internal_file: Internal path for file inside ZIP archive
    :param the_external_file: Full path and filename for file to copy to.
    :return: nothing
    """
    zip_file: zipfile.ZipFile
    read_zip_file: io.TextIOWrapper
    with zipfile.ZipFile(the_zip_path, 'r') as zip_file:
        if the_internal_file in zip_file.namelist():
            with zip_file.open(the_internal_file, mode="r") as read_zip_file:
                with open(the_external_file, 'wb') as out_file:
                    shutil.copyfileobj(read_zip_file, out_file)


def convert_int_to_str(the_int: int) -> str:
    """
    return string for int -- unless int is less than zero, then return empty string
    :param the_int: int to convert
    :return: return string for int -- unless int is less than zero, then return empty string
    """
    result: str = "NaN"
    if the_int >= 0:
        result = str(the_int)
    return result


def convert_str_to_int(the_str: str) -> int:
    """
    return int for string -- unless int is less than zero, then return empty string
    :param the_str: str to convert
    :return: return string for int -- unless int is less than zero, then return empty string
    """
    result: int = math.nan
    if "NaN" != the_str:
        result = int(the_str)
    return result


def convert_str_to_float(the_str: str) -> float:
    """
    return int for string -- unless int is less than zero, then return empty string
    :param the_str: str to convert
    :return: return string for int -- unless int is less than zero, then return empty string
    """
    result: float = math.nan
    if "NaN" != the_str:
        result = float(the_str)
    return result


def index_string_to_dict_str_str(the_str: str) -> dict[str, str]:
    """
    split index string key value list into dictionary.
    May have double-quotes around values.
    Values are strings.
    :param the_str: comma delimited key:value list
    :return: dictionary of key/values
    """
    my_dict: dict[str, str] = {}
    if '' != the_str:
        key_val_list: List[str] = the_str.strip().split(',')
        key_val: str
        for key_val in key_val_list:
            key_val = key_val.strip('"')
            splitted: List[str] = key_val.strip().split(':')
            key: str = splitted[0].strip()
            val: str = splitted[1].strip()
            my_dict[key] = val
    return my_dict


def index_string_to_dict_str_int(the_str: str) -> dict[str, int]:
    """
    split index string key value list into dictionary.
    May have double-quotes around values.
    Values are integers.
    :param the_str: comma delimited key:value list
    :return: dictionary of key/values
    """
    my_dict: dict[str, int] = {}
    if '' != the_str:
        key_val_list: List[str] = the_str.strip().split(',')
        key_val: str
        for key_val in key_val_list:
            key_val = key_val.strip('"')
            splitted: List[str] = key_val.strip().split(':')
            key: str = splitted[0].strip()
            val: str = splitted[1].strip()
            my_dict[key] = convert_str_to_int(val)
    return my_dict


def index_string_to_dict_str_float(the_str: str) -> dict[str, float]:
    """
    split index string key value list into dictionary.
    May have double-quotes around values.
    Values are Strings.
    :param the_str: comma delimited key:value list
    :return: dictionary of key/values
    """
    my_dict: dict[str, float] = {}
    if '' != the_str:
        key_val_list: List[str] = the_str.strip().split(',')
        key_val: str
        for key_val in key_val_list:
            key_val = key_val.strip('"')
            splitted: List[str] = key_val.strip().split(':')
            key: str = splitted[0].strip()
            val: str = splitted[1].strip()
            my_dict[key] = convert_str_to_float(val)
    return my_dict


def remove_unsafe_characters(the_str: str) -> str:
    """
    Make sure string is actually a string - Excel code started returning non-strings.
    Remove tabs, newlines, and line feeds from string.
    :param the_str: string to clean (might not even be a string)
    :return: cleaned string
    """
    # Excel processing sometimes
    the_str = str(the_str)
    # remove tab, line feeds and returns and double quotes
    the_str = re.sub("[\t\n\r\"]", "", the_str)
    return the_str
