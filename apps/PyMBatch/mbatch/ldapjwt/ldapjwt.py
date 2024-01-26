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

from typing import Dict, Optional
import json
import cryptography.fernet
import urllib3
import requests


def generate_key() -> str:
    key_b: bytes = cryptography.fernet.Fernet.generate_key()
    key_s: str = key_b.decode('utf-8')
    return key_s


def encrypt_string(the_str: str, the_key: str) -> str:
    key_b: bytes = the_key.encode('utf-8')
    str_b: bytes = the_str.encode('utf-8')
    fernet_obj: cryptography.fernet.Fernet = cryptography.fernet.Fernet(key_b)
    encrypted_b: bytes = fernet_obj.encrypt(str_b)
    encrypted_s: str = encrypted_b.decode('utf-8')
    return encrypted_s


def decrypt_string(the_str: str, the_key: str) -> str:
    key_b: bytes = the_key.encode('utf-8')
    str_b: bytes = the_str.encode('utf-8')
    fernet_obj: cryptography.fernet.Fernet = cryptography.fernet.Fernet(key_b)
    decrypted_b: bytes = fernet_obj.decrypt(str_b)
    decrypted_s: str = decrypted_b.decode('utf-8')
    return decrypted_s


def login_get_token(the_url_base: str, the_username: str, the_password: str) -> Optional[str]:
    """
    Authenticate and get a token for the user
    :param the_url_base: URL for API for LDAP-JWT
    :param the_username: username as a string
    :param the_password: password as a string
    :return: token as a string
    """
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    auth_url: str = f"{the_url_base}/authenticate"
    print(f"auth_url={auth_url}", sep='\n', flush=True)
    auth_data: Dict[str, str] = {"username": the_username, "password": the_password}
    # timeout is in seconds - wail up to 10 minutes (in case of reboot)
    auth_response: requests.Response = requests.post(auth_url, data=auth_data, verify=False, timeout=600)
    print(f"auth status code: {auth_response.status_code}", sep='\n', flush=True)
    result: str = ""
    if auth_response.status_code == 200:
        token: str = auth_response.json()["token"]
        print("auth response good", sep='\n', flush=True)
        print("auth response headers:", sep='\n', flush=True)
        print("Hit verify endpoint using new token", sep='\n', flush=True)
        verify_url: str = f"{the_url_base}/verify"
        print(f"verify_url={verify_url} ", sep='\n', flush=True)
        verify_data: Dict[str, str] = {"token": token}
        # timeout is in seconds - wail up to 10 minutes (in case of reboot)
        verify_response: requests.Response = requests.post(verify_url, data=verify_data, verify=False, timeout=600)
        print(f"verify status code: {verify_response.status_code}", sep='\n', flush=True)
        if verify_response.status_code == 200:
            print(f"verify response text: {json.dumps(verify_response.json(), indent=10, sort_keys=True)}", sep='\n', flush=True)
            result = token
        else:
            print(f"verify failed: {json.dumps(verify_response.json(), indent=10, sort_keys=True)}", sep='\n', flush=True)
    else:
        print("Bad response from auth endpoint")
    return result
