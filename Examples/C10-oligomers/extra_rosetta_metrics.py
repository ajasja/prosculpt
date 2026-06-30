import os
import subprocess
import hashlib
import sys
# from IPython.utils import io
# import py3Dmol
import pyrosetta

import pandas as pd

basedir = '/home/aljubetic/projects/2024-03-15__oligomer_selections/fastas'
df = pd.read_csv('/home/akonstantinova/shared/18_03_2024_oligomers/best_with_delta_lddt_9_11.csv')
ids = list(df['ID'])
print(ids)

pyrosetta.init("-ignore_unrecognized_res -ignore_zero_occupancy  -corrections::beta_nov16 true")

protocol = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects().create_from_string(
"""
<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="sfxn_clean" weights="beta_nov16" />
    </SCOREFXNS>

    <RESIDUE_SELECTORS>
        <Chain name="chainA" chains="A" />
        <Chain name="chainB" chains="B" />
        <ScoreTermValueBased name="clashing_res" score_type="fa_rep" score_fxn="sfxn"  lower_threshold="3" upper_threshold="99999" />


    </RESIDUE_SELECTORS>
    <TASKOPERATIONS>
        <IncludeCurrent name="ic" /> //includes input pdbs rotamers
        <LimitAromaChi2 name="limitaro" chi2max="110" chi2min="70" include_trp="1" /> //disallow extreme aromatic rotamers
        <ExtraRotamersGeneric name="ex1_ex2" ex1="1" ex2="1" /> //use ex1 ex2 rotamers
        <RestrictToRepacking name="repack_only" />  //for minimize/repack
    </TASKOPERATIONS>

    <FILTERS>
        <NetCharge name="chargeA" chain="1" confidence="0" />
        <ShapeComplementarity name="sc2" min_sc="0.6" verbose="1" quick="0" residue_selector1="chainA" residue_selector2="chainB" write_int_area="1" write_median_dist="1" confidence="0" />
        <ExposedHydrophobics name="exposed_hydrop" sasa_cutoff="20" threshold="0" confidence="0"/>
        <ContactMolecularSurface name="cms" distance_weight="0.5" use_rosetta_radii="true" apolar_target="0"
        target_selector="chainA" binder_selector="chainB" confidence="0" />

        <BuriedUnsatHbonds name="vbuns"  report_all_heavy_atom_unsats="true" scorefxn="sfxn_clean" ignore_surface_res="false" print_out_info_to_pdb="true" atomic_depth_selection="5" burial_cutoff="1000" residue_surface_cutoff="42.5" dalphaball_sasa="0" confidence="0"  only_interface="false"  />
	    	<BuriedUnsatHbonds name="sbuns" report_all_heavy_atom_unsats="true" scorefxn="sfxn_clean" cutoff="4" residue_surface_cutoff="42.5" ignore_surface_res="false" print_out_info_to_pdb="true" dalphaball_sasa="0" probe_radius="1.1" atomic_depth_selection="5.5" atomic_depth_deeper_than="false" only_interface="false" confidence="0" />
    </FILTERS>

    <MOVERS>
        <SwitchChainOrder name="take_two_chains" chain_order="12"/> #Or whatever chain you need to get one neighbor
        <MinMover name="minimize_sc_all" scorefxn="sfxn_clean" bb="0" chi="1" />
        <InterfaceAnalyzerMover name="analyze_interface" scorefxn="sfxn_clean"
        packstat="1" interface_sc="1" use_jobname="1"
        jump="1" scorefile_reporting_prefix="IA" />
        <DumpPdb name="dump_pose" fname="dump.pdb" />
        <ddG name="ddG_no_repack" translate_by="1000" scorefxn="sfxn_clean"  task_operations="repack_only,ic,ex1_ex2" relax_mover="minimize_sc_all"
            repack_bound="0"
            relax_bound="0"
            repack_unbound="0"
            relax_unbound="1"
        jump="1"
        dump_pdbs="0"   />
    </MOVERS>
    <SIMPLE_METRICS>
        <SapScoreMetric name="sap_score" />
        <SelectedResiduesPyMOLMetric name="clashing_res" residue_selector="clashing_res" custom_type="clashing_res" />
    </SIMPLE_METRICS>
    <PROTOCOLS>
        <Add mover="take_two_chains" />
        Add mover='dump_pose' />
        <Add mover="minimize_sc_all" />
        Add mover="analyze_interface" />
        <Add mover="ddG_no_repack" />
        <Add filter="chargeA" />
        <Add filter="exposed_hydrop" />
        <Add filter="cms" />
        <Add filter="vbuns" />
        <Add filter="sbuns" />
        <Add metrics="sap_score" />
        <Add metrics="clashing_res" />

        <Add filter="sc2" />
    </PROTOCOLS>
</ROSETTASCRIPTS>

""").get_mover("ParsedProtocol")

st = {}
for id in ids:
    try:
        for pdb in os.listdir(basedir+'/'+id):
            if 'relaxed' in pdb:
                print(pdb)
                pose = pyrosetta.pose_from_pdb(basedir+'/'+id+'/'+pdb)
                protocol.apply(pose)
                
                extra_dict = dict(pose.scores)
                job_pairs = pyrosetta.rosetta.protocols.jd2.get_string_real_pairs_from_current_job() 

                for name in job_pairs:
                    if  not name in extra_dict: 
                        extra_dict[name] = job_pairs[name]

                print(extra_dict)
                st[pdb]=extra_dict
    
    except FileNotFoundError:
        st[id]={}

import json
with open('rosettastats_9_10_11_rstats.json', 'w') as f:
    json.dump(st, f)