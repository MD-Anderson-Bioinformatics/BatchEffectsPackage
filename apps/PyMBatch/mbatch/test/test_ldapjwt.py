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


import unittest
import os
import mbatch.ldapjwt.ldapjwt
import mbatch.test.common


# pylint: disable=too-many-instance-attributes
class TestJob(unittest.TestCase):
    """
    Class for setting up Job testing - clear/make directory for output
    """

    # do not set method variables, as they should be initialized in the init function
    # No local method variables

    def setUp(self: 'TestJob') -> None:
        """
        setup script to clear and re-populate test directory
        :return:
        """
        #############################
        # files for configout testing
        # if self._testMethodName == 'test_process_configout_dir':
        #     print(f"TestJob::setUp dynamic_test_dir={dynamic_test_job_dir}", flush=True)
        #     if os.path.exists(dynamic_test_job_dir):
        #         shutil.rmtree(dynamic_test_job_dir)
        #     os.makedirs(dynamic_test_job_dir)
        #     print(f"TestJob::setUp static_test_dir={static_test_job_dir}", flush=True)
        #############################
        print("TestJob:setUp done", flush=True)

    def test_crypt(self: 'TestJob') -> None:
        """
        test the create_index_archive function with the configout test directory
        :return: nothing
        """
        print("test_crypt", flush=True)
        test_key: str = mbatch.ldapjwt.ldapjwt.generate_key()
        print(f"test_key={test_key}", sep='\n', flush=True)
        original_str: str = "ABCDEFG01234567890~!@#$"
        print(f"original_str={original_str}", sep='\n', flush=True)
        encrypted_str: str = mbatch.ldapjwt.ldapjwt.encrypt_string(original_str, test_key)
        print(f"encrypted_str={encrypted_str}", sep='\n', flush=True)
        decrypted_str: str = mbatch.ldapjwt.ldapjwt.decrypt_string(encrypted_str, test_key)
        print(f"decrypted_str={decrypted_str}", sep='\n', flush=True)
        assert original_str == decrypted_str, "original and decrypted strings are different"

    def test_ldap_jwt(self: 'TestJob') -> None:
        # these files are only readable by application's development user
        if os.path.exists('/BEA/key.txt'):
            test_key: str = mbatch.test.common.read_file_to_string('/BEA/key.txt')
            encrypted_uid: str = mbatch.test.common.read_file_to_string('/BEA/uid.txt')
            encrypted_pwd: str = mbatch.test.common.read_file_to_string('/BEA/pwd.txt')
            decrypted_uid: str = mbatch.ldapjwt.ldapjwt.decrypt_string(encrypted_uid, test_key)
            decrypted_pwd: str = mbatch.ldapjwt.ldapjwt.decrypt_string(encrypted_pwd, test_key)
            ldap_jwt_url: str = mbatch.test.common.read_file_to_string('/BEA/ldapjwt.txt')
            token: str = mbatch.ldapjwt.ldapjwt.login_get_token(ldap_jwt_url, decrypted_uid, decrypted_pwd)
            print(f"{token}", sep='\n', flush=True)
            assert '' != token, "token not retrieved"
        else:
            print("test_ldap_jwt not run, no test files", flush=True)

# pylint: enable=too-many-instance-attributes


if __name__ == '__main__':
    unittest.main()
