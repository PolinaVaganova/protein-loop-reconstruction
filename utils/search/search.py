# flake8: noqa

import json
import pathlib
import re
import shlex
import subprocess
from os.path import isfile

import gemmi
import requests
from Bio.PDB import PDBList
from tqdm import tqdm


def get_column_label(mtz):
    pattern_str = r"[IF]-?(OBS)?(\([-+]\))?"
    pattern = re.compile(pattern_str, re.IGNORECASE)
    for label in mtz.column_labels():  # type: gemmi.Mtz.Column
        if pattern.fullmatch(label):
            return label


class ProcessPDB:
    def __init__(self):
        self.pdb_codes = []
        self.cif_list = []
        self.no_altlocs_pdbs = []
        self.pdb_codes_no_twinning = []
        # If you wish to run only the last step to filter more by resolution. In this case leave only self.prep_mtzs() in run()
        # self.pdb_codes_no_twinning = ['2qdo', '5lhx', '6pcf', '5teo', '4o8b', '6rlr', '3shu', '5xbh', '6kmq', '7mpw',
        #                               '5zcx', '2nn4', '2o89', '2o87', '5f6a', '2o85', '5tut', '7ccd', '1kem', '1smt',
        #                               '2eql', '1ail', '1k40', '5hu4', '3le4', '4e5r', '4r22', '3cro', '4x37', '3k9p',
        #                               '2fht', '4zh9', '1x6j', '1q2y', '4f26', '1anu', '2aak', '7kl2', '1dt4', '2ont',
        #                               '1z27', '6rjx', '4wfw', '6a2h', '3rd3', '6cyr', '3dvp', '6fc5', '7mbk', '1du5',
        #                               '2snw', '4niq', '1ov9', '3hpm', '3ndf', '4fis', '1ytg', '3fis', '149l', '5x56',
        #                               '5h79', '2i0x', '1thu', '2dff', '6njg', '3jsv', '3u71', '5k9n', '5en1', '1g83',
        #                               '5uh7', '2auj', '1ae6', '5bn0', '2rbf', '6i9j', '4e9h', '1i1f', '1eez', '1i7t',
        #                               '6ekm', '7df7', '7bzd', '1lou', '3da5', '4rlm', '5bw5', '4hv2', '1z91', '5oq0',
        #                               '7c20', '4rln', '2d1v', '2epe', '4e9g', '5ju5', '6vjw', '4i48', '5ego', '5xzx',
        #                               '6gf0', '2yvh', '7cnb', '3ri4', '6ul3', '1uef', '1wmi', '1jhe', '5i5o', '5y2e',
        #                               '2e8c', '4yyg', '1jhc', '1ab0', '7lke', '4zm2', '6ebb', '1a1e', '1a1c', '5ifn',
        #                               '5ep5', '1izz', '1a07', '1a08', '2qc7', '6agw', '5l0o', '1bt6', '3q5x', '2h6c',
        #                               '4zyd', '4qmf', '5eg0', '4j6q', '5nod', '1xnh', '5yah', '4ozk', '4opa', '5eep',
        #                               '4lh9', '6pse', '6p7t', '3bjv', '7e7g', '1wjg', '4hsv', '5jjz', '6t2e', '6j62',
        #                               '6qlz', '2pwj', '5yry', '4iyl', '1d0e', '4pdb', '4rsv', '6xwu', '1bu2', '5lx2',
        #                               '2h1g', '2q2p', '1ayb', '6jb6', '1squ', '4noj', '3ong', '4ka4', '6n5m', '1ao6',
        #                               '4pu8', '5dpr', '7dy9', '3sgo', '6l2b', '4mo0', '7d35', '3ofj', '6j8f', '6ls6',
        #                               '6az5', '5lux', '2b8i', '4fdx', '167l', '5ejo', '6gbq', '3fk9', '2h1k', '5yac',
        #                               '3nth', '3toe', '4qxz', '3usv', '5by2', '1kxc', '1czg', '4qr1', '1kxe', '1kxa',
        #                               '1kxd', '1kxf', '1kxb', '5u9k', '6sbw', '4f4m', '6m5i', '4u7y', '5j4g', '1zw3',
        #                               '3g2w', '5yv7', '1xwj', '3r3p', '5u9a', '1xhm', '3kxe', '6hnk', '3lz2', '1c2a',
        #                               '4et7', '5eri', '1oz9', '1xwc', '5efx', '4qr0', '4wbf', '4om7', '2fxk', '4euw',
        #                               '5yb4', '2huz', '4tr5', '1pdb', '7m5u', '2fpf', '1ll6', '1ar2', '5lty', '1v1p',
        #                               '5yzj', '6lmj', '1jbg', '5yzl', '2ldx', '2fzf', '2z0w', '2x5s', '4uf2', '5a0g',
        #                               '7myf', '1hij', '1hik', '3ggr', '6jhe', '6yqf', '6l2s', '4e6s', '4uf3', '1mqa',
        #                               '4m8q', '4r8e', '2if4', '6l9u', '4odx', '1kct', '6ign', '6w69', '5w35', '5f22',
        #                               '6w67', '4ezq', '4lon', '1zr5', '6wiq', '4lp2', '1c3g', '1a3f', '5f5j', '3hnr',
        #                               '177l', '5f5g', '6lfb', '7kkf', '6a2b', '5uic', '4lmz', '2d2q', '1aqt', '5z2t',
        #                               '4z8m', '4is6', '4iuf', '1hkl', '2no2', '7d3t', '1mam', '5oc6', '3avz', '3hw2',
        #                               '2rae', '4s2q', '3vyy', '5qkh', '1pgz', '5qlp', '4kz1', '5bxf', '1qqd', '1ay9',
        #                               '6h97', '1tc3', '5wwi', '5wxd', '5yb9', '1y6e', '3gtn', '1tfo', '5xmm', '2qyp',
        #                               '4nod', '1u0o', '5gtx', '6xh1', '3lb9', '1bz9', '1jpg', '7da2', '5xdm', '4odd',
        #                               '7aye', '1k7a', '2f15', '5wwu', '6jq4', '5k1y', '2np2', '3voe', '4mji', '4lgs',
        #                               '3vod', '6f5b', '6fsp', '3q0x', '6g1l', '6cf2', '6jh1', '5egk', '7c10', '5znd',
        #                               '6f8f', '5y7i', '1q1b', '1g6y', '6gp9', '6e8c', '4nqs', '3g9z', '3hfi', '5cn2',
        #                               '1osz', '5exv', '4yg4', '4i6p', '2lve', '1wvi', '3ilh', '4f2h', '4r0b', '5ev7',
        #                               '6m2k', '1frg', '3o1b', '3h1b', '4l2w', '6gy3', '3ee8', '3x0w', '1ttw', '1a7h',
        #                               '6ulr', '7njh', '1fp9', '4zou', '3thf', '7nji', '1t3b', '6m4c', '6euj', '6y93',
        #                               '3gjq', '1mdm', '3e4h', '7e90', '6iju', '2irt', '7pq6', '5ic6', '5hm1', '4k8u',
        #                               '4om2', '7jjl', '3jrf', '4n2s', '5e3n', '4nsr', '3iv5', '3jri', '7cyg', '3jrg',
        #                               '5dtd', '3jrd', '3jrh', '5e3l', '6vsk', '1itb', '5ds9', '3s4u', '5e3o', '1uhd',
        #                               '5mjw', '2zkf', '5e3m', '3jrc', '6hq2', '3m7m', '3jra', '3jrb', '3jre', '3rnv',
        #                               '1qz7', '4n7v', '5dhx', '4ob4', '1kql', '5o1z', '6udg', '4db4', '3va2', '4kji',
        #                               '5zz3', '7vjm', '1pdr', '4n2q', '1ibc', '3mfk', '2dbo', '5toi', '1ry7', '2nue',
        #                               '6vmt', '2zke', '5no6', '1hyq', '7rg4', '4i6z', '4wxm', '6jsx', '3btr', '4fq0',
        #                               '1f36', '4p39', '7dhb', '1vll', '7dh4', '6yl2', '4h22', '7dwo', '4htv', '6ti2',
        #                               '3s1b', '1ri7', '3lg8', '4lg2', '3ve6', '3ojb', '3m8v', '4xqq', '4nkj', '5x8y',
        #                               '3s63', '2zte', '4mvb', '5ekg', '4enj', '2haf', '1py1', '7dgx', '6yek', '1ooc',
        #                               '1zkz', '4enm', '6i5e', '4nbo', '5i4t', '1z3g', '1sv6', '5zvm', '5xuo', '1tgk',
        #                               '3f6l', '4gjh', '4eq2', '3d4v', '2etn', '4yg1', '6je7', '6nrl', '6pln', '5hlt',
        #                               '2owy', '4pts', '1agx', '3cwu', '3w6v', '5hwf', '1lo5', '1bwz', '5kk1', '4l67',
        #                               '3cvt', '6wk5', '7ctn', '6gef', '6ri3', '3cw7', '4j2l', '6m10', '6l2r', '6ryi',
        #                               '6kbm', '3m7g', '1isn', '4iht', '2emt', '5h77', '1bcm', '4trz', '4i1t', '4pba',
        #                               '6tpt', '3mtv', '4try', '2h1o', '1f9j', '6wc5', '1a9n', '7byj', '4trw', '4y0f',
        #                               '3hon', '3gm1', '4l22', '5inc', '6w6f', '2yw7', '5m2e', '5ygj', '6h02', '6s9s',
        #                               '1czd', '2io5', '1edz', '1b8h', '2z8h', '4rwq', '3gpn', '3te3', '1imh', '1a9b',
        #                               '6eo2', '5hpc', '3c0l', '4k48', '6n3i', '3rqg', '2f8s', '1zyl', '7b22', '1fsk',
        #                               '6w4b', '7eu2', '6h9n', '5vfx', '3wzf', '209l', '2cwe', '3hfm', '2z3x', '3btz',
        #                               '2ems', '3gm2', '1kq2', '4dm9', '3bis', '4is7', '5l1m', '2gjw', '6sbx', '2qdq',
        #                               '6t3h', '2i9t', '6o03', '4zus', '2f54', '5x9b', '5x07', '7f4w', '2pjy', '3h2t',
        #                               '3vzr', '4zuw', '2ge8', '1iaw', '5k23', '6wus', '1t2k', '6s9r', '6jij', '3nyl',
        #                               '4xrk', '1b6u', '4gtu', '3gu0', '3q5f', '3wzi', '6jzd', '1owr', '1otp', '4wuu',
        #                               '5k22', '3omw', '6smd', '5vxk', '5v6h', '4oph', '2nrs', '4fgn', '6eba', '6alx',
        #                               '5vu6', '4enn', '6ajk', '6f5f', '4u8g', '5m8a', '5yzz', '4zxd', '2z0s', '5uq2',
        #                               '2fse', '5m8j', '6j4r', '4ih4', '5m8k', '3qjp', '4dg7', '2ihr', '5hmo', '6ov0',
        #                               '1xfb', '6mqm', '6vg2', '4dt0', '7ccj', '1x8z', '3txq', '1tij', '6yxj', '6yhl',
        #                               '1gc7', '6uk4', '7e74', '5dbk', '5c3n', '1c9i', '1c9l', '4r9y', '6c2s', '6kfp',
        #                               '6jh0', '4dzs', '6jn2', '4j6s', '5bz2', '6uk2', '6gc9', '4eoz', '5zmc', '6l9w',
        #                               '3oqn', '6vj9', '5kh1', '5xgl', '4l5i', '3bba', '3tkn', '5xby', '5xoc', '7bme',
        #                               '3uj3', '5o7g', '7lio', '3ma5', '5dv0', '3hzr', '2rh5', '7n94', '1tx9', '6jw3',
        #                               '1d5f', '6jw4', '1ryx', '6jw2', '6jw5', '7e16', '4qqb', '3qjj', '4y21', '2zr3',
        #                               '5m9e', '1p4l', '1f02', '6o1z', '5mu7', '3kli', '5u3j', '1td6', '5b61', '6ysh',
        #                               '4fi9', '4kr9', '4tzv', '3vrl', '3r3h', '6kza', '6ozz', '6ozx', '4fq7', '7bz6',
        #                               '6lfp', '1n9r', '4d9v', '1qb3', '4mbr', '4qpq', '4qx2', '4qx0', '2qyg', '1tfc',
        #                               '5d8w', '1kv6', '4qx3', '4qx1', '7bpk', '6oz8', '1tzi', '6lvp', '4rbo', '3ks2',
        #                               '4mjn', '5ujy', '4zw0', '2grl', '4tui', '6kcr', '5jp3', '3a7p', '4mt7', '3wtp',
        #                               '6o0s', '6sws', '3oh9', '6nly', '2gqq', '6amk', '5djn', '2yvc', '2io2', '7cjn',
        #                               '7d9m', '2io3', '3vgt', '3do7', '7f9n', '3ogd', '7cd6', '3oh6', '5uco', '7f9m',
        #                               '1ow7', '1ow8', '5x47', '1qso', '2qrx', '1k05', '6o6v', '1dov', '3hgk', '6o71',
        #                               '5jly', '1ytd', '6o6s', '5bsa', '6amn', '4e79', '1vbx', '1c1g', '2fkd', '1u78',
        #                               '1sj4', '1ulq', '3n8e', '1fcc', '2dpd', '5ay9', '6oq6', '5gtb', '1xni', '3qjl',
        #                               '4ydz', '104l', '2av5', '3rb9', '1zr4', '4e5z', '1av1', '3prd', '4j40', '3kg5',
        #                               '5c76', '2aw6', '6igv', '4kw3', '3b3l', '7cd5', '5x6d', '6vfm', '3mfm', '3mqk',
        #                               '6tyd', '5hyp', '2jgz', '6h41', '6ajl', '6ef4', '2f86', '4ixp', '6aj4', '6hs6',
        #                               '3ooc', '1u2j', '1i1g', '6zne', '7dc6', '6mqs', '3qoq', '1f7v', '3oq9', '4rzi',
        #                               '6q5v', '1n9s', '4q0c', '6iak', '3qtn', '1cax', '6hb3', '5lfj', '5xm9', '4zyh',
        #                               '5zli', '3tnu', '5ila', '3a2k', '1zcf', '2hxy', '2zl3', '1mc4', '2dwm', '2zl2',
        #                               '2qvj', '4twh', '4u4a', '4iq4', '6vg8', '5tbk', '7dtg', '5ept', '5z3r', '3hhg',
        #                               '6x1i', '1wmk', '6j7g', '1x9g', '1cii', '1otc', '6l93', '4wl2', '1pci', '6ama',
        #                               '5xhv', '5dy0', '2o6g', '4dej', '4nn1', '1krq', '2pek', '2pem', '2pen', '2pej',
        #                               '2z46', '4zsv', '3ibb', '1ngm', '6xwt', '4hrn', '7cqj', '5n47', '2yyb', '6wl2',
        #                               '3hhz', '6wl3', '1koh', '3jvz', '3jw0', '6wl4', '6x0r', '1koo', '3aae', '1z1b',
        #                               '3kml', '3agr', '4qyj', '5ctr', '3s5c', '4m6f', '3t1p', '3mkr', '6r2s', '5x3f',
        #                               '5j7t', '3oed', '6rmm', '4wsc', '5oxv', '4zjd', '3uw0', '2hxb', '2nwc', '5opw',
        #                               '3f9v', '5cdi', '7asc', '4ijs', '6wg6', '7fbd', '5ch6', '4qrm', '6zn9', '7d8t',
        #                               '4omt', '7bwk', '6c6d', '6ki7', '4eld', '5eeb', '4d9j', '4itv', '3h4p', '6job',
        #                               '6j4a', '2du7', '4ro0', '1pma', '6ahq', '5v5f', '5yk3', '3dyu', '5nl0', '5jhc',
        #                               '5xrc', '5l5n', '6fdb', '5mnt', '3a6m', '4zln', '4zlj', '4zll', '4qvg', '4wgl',
        #                               '1oqe', '1oqd', '1jh5', '4ypi', '6m3c', '2iir', '3wof', '3win', '4i9c', '1eaa',
        #                               '1dpc', '1dpb', '3p0s', '6okq', '4k7h', '5lhy', '3rqc', '6hkt', '5l5m', '2eu1',
        #                               '6ria', '1gav', '5c73', '1m06', '1a6c', '4u6u', '1dnv', '2fz1', '4zpy', '2g34',
        #                               '2qij', '2ft1']

        self.final_list = []
        self.exceptions_cif_as_mtz = [
            # '1xnh', # this one needed manual processing with --incompatible_flags_to_work_set flag
            # '2no2', # this one needed manual processing with --symmetry "P 42 21 2" flag
        ]
        self.session = requests.Session()

    def run(self):
        self.search4pdbs()
        self.download_pdbs()
        self.filter_small_structures()
        self.no_altloc_structures()
        self.download_reports()
        self.check_reports()
        self.to_refine()
        self.prep_mtzs()

    def search4pdbs(self):
        """
        condition 1: contains proteins, NAs or their hybrids
        condition 2: has experimental data and x-ray
        condition 3: no ligands, no more than 50 water molecules, no modified residues
        condition 4: no membrane proteins by keywords
        condition 4: no membrane proteins by annotations
        :return:
        """
        url = "https://search.rcsb.org/rcsbsearch/v1/query"

        json_request = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "group",
                        "nodes": [
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                                    "operator": "exact_match",
                                    "negation": False,
                                    "value": "Nucleic acid (only)",
                                },
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                                    "operator": "exact_match",
                                    "negation": False,
                                    "value": "Protein (only)",
                                },
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                                    "operator": "exact_match",
                                    "negation": False,
                                    "value": "Protein/NA",
                                },
                            },
                        ],
                        "logical_operator": "or",
                    },
                    {
                        "type": "group",
                        "nodes": [
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "exptl.method",
                                    "operator": "exact_match",
                                    "negation": False,
                                    "value": "X-RAY DIFFRACTION",
                                },
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "rcsb_accession_info.has_released_experimental_data",
                                    "operator": "exact_match",
                                    "negation": False,
                                    "value": "Y",
                                },
                            },
                        ],
                        "logical_operator": "and",
                    },
                    {
                        "type": "group",
                        "nodes": [
                            {
                                "type": "group",
                                "logical_operator": "and",
                                "nodes": [
                                    {
                                        "type": "terminal",
                                        "service": "text",
                                        "parameters": {
                                            "attribute": "rcsb_polymer_entity_feature_summary.count",
                                            "operator": "equals",
                                            "negation": False,
                                            "value": 0,
                                        },
                                    },
                                    {
                                        "type": "terminal",
                                        "service": "text",
                                        "parameters": {
                                            "attribute": "rcsb_polymer_entity_feature_summary.type",
                                            "operator": "exact_match",
                                            "value": "modified_monomer",
                                            "negation": False,
                                        },
                                    },
                                ],
                                "label": "nested-attribute",
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "rcsb_entry_info.deposited_solvent_atom_count",
                                    "operator": "less_or_equal",
                                    "negation": False,
                                    "value": 50,
                                },
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "rcsb_entry_info.deposited_nonpolymer_entity_instance_count",
                                    "operator": "equals",
                                    "negation": False,
                                    "value": 0,
                                },
                            },
                            # {
                            #     "type": "terminal",
                            #     "service": "text",
                            #     "parameters": {
                            #         "attribute": "entity_poly.rcsb_mutation_count",
                            #         "operator": "equals",
                            #         "negation": False,
                            #         "value": 0
                            #     }
                            # }
                        ],
                        "logical_operator": "and",
                    },
                    {
                        "type": "group",
                        "nodes": [
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "struct_keywords.text",
                                    "operator": "contains_words",
                                    "negation": True,
                                    "value": "membrane",
                                },
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "struct_keywords.pdbx_keywords",
                                    "operator": "contains_words",
                                    "negation": True,
                                    "value": "membrane",
                                },
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "rcsb_uniprot_protein.name.value",
                                    "operator": "contains_words",
                                    "negation": True,
                                    "value": "membrane",
                                },
                            },
                        ],
                        "logical_operator": "and",
                    },
                    {
                        "type": "group",
                        "logical_operator": "and",
                        "nodes": [
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "negation": True,
                                    "attribute": "rcsb_polymer_entity_annotation.type",
                                    "operator": "exact_match",
                                    "value": "PDBTM",
                                },
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "negation": True,
                                    "attribute": "rcsb_polymer_entity_annotation.type",
                                    "operator": "exact_match",
                                    "value": "MemProtMD",
                                },
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "negation": True,
                                    "attribute": "rcsb_polymer_entity_annotation.type",
                                    "operator": "exact_match",
                                    "value": "OPM",
                                },
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "negation": True,
                                    "attribute": "rcsb_polymer_entity_annotation.type",
                                    "operator": "exact_match",
                                    "value": "mpstruc",
                                },
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "negation": True,
                                    "attribute": "rcsb_polymer_entity_annotation.type",
                                    "operator": "exact_match",
                                    "value": "mpstruc",
                                },
                            },
                        ],
                    },
                ],
                "label": "text",
            },
            "return_type": "entry",
            "request_options": {
                "pager": {
                    "start": 0,
                    "rows": 6000,  # this should be enough with overhead
                },
                "scoring_strategy": "combined",
                "sort": [
                    {
                        "sort_by": "rcsb_entry_container_identifiers.entry_id",
                        "direction": "asc",
                    }
                ],
            },
        }
        print("querying PDB...\n")
        report = self.session.get(f"{url}?json={json.dumps(json_request)}")

        if report.status_code != 200:
            print(f"return code: {report.status_code}")
            return False

        result = report.json()
        if result:
            print("Found number of PDB entries:", result["total_count"])
            self.pdb_codes = sorted([r["identifier"] for r in result["result_set"]])
        else:
            print("Failed to retrieve results")
            self.pdb_codes = []

    def download_pdbs(self, dummy_run=True):
        pdbdownload = PDBList(verbose=False)
        if dummy_run:
            for code in tqdm(self.pdb_codes, desc="Download PDBs"):
                if not isfile(f"db/pdb{code.lower()}.ent"):
                    self.cif_list.append(code)
            return

        for code in tqdm(self.pdb_codes, desc="Download PDBs"):
            pdbdownload.retrieve_pdb_file(
                code,
                obsolete=False,
                pdir="db/",
                file_format="pdb",
                overwrite=False,
            )
            if not isfile(f"db/pdb{code.lower()}.ent"):
                self.cif_list.append(code)
                pdbdownload.retrieve_pdb_file(
                    code,
                    obsolete=False,
                    pdir="db/",
                    file_format="mmCif",
                    overwrite=False,
                )

    def filter_small_structures(self, min_side_length=16):
        volumes = {}
        for code in tqdm(self.pdb_codes, desc="Cut off by size"):
            if code in self.cif_list:
                pdb = gemmi.read_structure(f"db/{code.lower()}.cif")
            else:
                pdb = gemmi.read_pdb(f"db/pdb{code.lower()}.ent")
            if (
                pdb.cell.a > min_side_length
                and pdb.cell.b > min_side_length
                and pdb.cell.c > min_side_length
            ):  # adjust these values
                volumes[code] = pdb.cell.volume
        sorted_by_value = sorted(volumes.items(), key=lambda kv: kv[1])
        self.not_so_small_pdbs = [k[0] for k in sorted_by_value]

    def check_occupancy_and_missing_residues(self, input_structure):
        input_structure.remove_waters()
        input_structure.remove_alternative_conformations()
        input_structure.remove_empty_chains()
        for model in input_structure:
            for chain in model:
                chain_range = []
                for residue in chain:
                    chain_range.append(residue.seqid.num)
                    for atom in residue:  # type: gemmi.Atom
                        if atom.occ < 1.0:
                            return False
                for i in range(min(chain_range), max(chain_range)):
                    if i not in chain_range:
                        return False
        return True

    def no_altloc_structures(self):
        for code in tqdm(
            self.not_so_small_pdbs, desc="Cut off altloc and gaps in chains"
        ):
            if code in self.cif_list:
                pdb = gemmi.read_structure(f"db/{code.lower()}.cif")
            else:
                pdb = gemmi.read_pdb(f"db/pdb{code.lower()}.ent")
            if self.check_occupancy_and_missing_residues(pdb):
                self.no_altlocs_pdbs.append(code)

    def download_report(self, pdb_code):
        filepath = f"db/pdb{pdb_code}_report.pdf"
        if isfile(filepath):
            return True
        link = f"http://files.rcsb.org/pub/pdb/validation_reports/{pdb_code[1:3]}/{pdb_code}/{pdb_code}_full_validation.pdf"
        report = self.session.get(link)
        if report.status_code != 200:
            return False
        with open(filepath, "wb") as f:
            f.write(report.content)
        return True

    def download_reports(self):
        for code in tqdm(self.no_altlocs_pdbs, desc="Download validation reports"):
            if not self.download_report(code.lower()):
                print(code, "couldn't download report")

    def check_report(self, pdb_code):
        s = subprocess.run(
            shlex.split(f"pdfgrep twin db/pdb{pdb_code}_report.pdf"),
            stdout=subprocess.PIPE,
        )
        report = s.stdout.decode("ascii")
        if "No twinning" in report:
            self.pdb_codes_no_twinning.append(pdb_code)

    def check_reports(self):
        for code in tqdm(self.no_altlocs_pdbs, desc="Check validation reports"):
            self.check_report(code.lower())
        # print('self.pdb_codes_no_twinning = ')
        # print(self.pdb_codes_no_twinning)

    def download_sfs(self, pdb_code):
        filepath = f"structure_factors/{pdb_code}-sf.cif"
        if isfile(filepath):
            return True
        link = f"https://files.rcsb.org/download/{pdb_code}-sf.cif"
        report = self.session.get(link)
        if report.status_code == 200:
            with open(filepath, "wb") as f:
                f.write(report.content)
            return True
        else:
            return False

    def to_refine(self):
        pathlib.Path("structure_factors").mkdir(parents=True, exist_ok=True)
        for code in tqdm(self.pdb_codes_no_twinning, desc="Download structure factors"):
            if not (self.download_sfs(code)):
                print(f"Problems downloading {code}")

    def prep_mtzs(self):
        for code in tqdm(
            self.pdb_codes_no_twinning, desc="Prepare <code>_P1.mtz files"
        ):
            if code in self.exceptions_cif_as_mtz:
                continue
            filename_base = f"structure_factors/{code}-sf"
            if not isfile(f"{filename_base}.mtz"):
                command_cis_as_mtz = (
                    f"phenix.cif_as_mtz {filename_base}.cif --ignore_bad_sigmas --merge "
                    f"--output_file_name={filename_base}.mtz"
                )
                if code == "1xnh":
                    command_cis_as_mtz += " --incompatible_flags_to_work_set"
                if code == "2no2":
                    command_cis_as_mtz += ' --symmetry "P 42 21 2"'

                s = subprocess.run(
                    shlex.split(command_cis_as_mtz),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                report = s.stdout.decode("ascii")
                # print(report)
                report = s.stderr.decode("ascii")
                # print(report)
            report = ""
            mtz = gemmi.read_mtz_file(f"{filename_base}.mtz")
            data_label = get_column_label(mtz)
            if (
                "+" in data_label
            ):  # handle the case when anomalous scattering is present
                mod_label = data_label.replace("+", "-")
                data_label = f"{data_label},{mod_label}"
            if not isfile(f"{filename_base}_P1.mtz"):
                command_line = f'phenix.reflection_file_converter --expand_to_p1 {filename_base}.mtz --write_mtz_amplitudes --mtz_root_label="FOBS" --label={data_label} --generate_r_free_flags --non_anomalous --mtz {filename_base}_P1.mtz'
                s = subprocess.run(
                    shlex.split(command_line),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                report = s.stdout.decode("ascii")
                # print(report)
                report = s.stderr.decode("ascii")
                # print(report)
            if "Unknown --label=" in report:
                print(report)
            elif mtz.resolution_high() > 1.5:
                self.final_list.append(code)
            # break
        print("Final list of structures:")
        print(self.final_list)


if __name__ == "__main__":
    db = ProcessPDB()
    db.run()
