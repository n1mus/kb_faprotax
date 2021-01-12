# -*- coding: utf-8 -*-
from configparser import ConfigParser
import os
import logging
import time
import unittest
import uuid
from shutil import rmtree, copytree
import json
from unittest.mock import patch
import pandas as pd
from urllib.parse import urlparse
import re
from pytest import raises

from kb_faprotax.kb_faprotaxImpl import kb_faprotax
from kb_faprotax.kb_faprotaxServer import MethodContext
from kb_faprotax.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace


from kb_faprotax.util.dprint import dprint, where_am_i
from kb_faprotax.util.varstash import Var
from kb_faprotax.util.kbase_obj import GenomeSet, Genome
from mock import * # mock business
from upa import * # upa library
import config




class Test(config.BaseTest):

    def test_GenomeSet_input(self):
        ret = config.get_serviceImpl().run_FAPROTAX(
            config.ctx, 
            {
                **self.ws,
                'input_upa': refseq,
            }
        )

    def test_dup_GenomeSet(self):
        ret = config.get_serviceImpl().run_FAPROTAX(
            config.ctx, 
            {
                **self.ws,
                'input_upa': refseq_dup,
            }
        )



