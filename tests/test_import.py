import pytest

import os
import sys

import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html

# Add the directory containing the module to sys.path
# root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin'))
# sys.path.append(root_dir)
root_dir = '/Users/danielbergman/PhysiCell-Studio_git/bin'
sys.path.append(root_dir)

# Now you can import the module
import studio
from pretty_print_xml import pretty_print


from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication

class StudioTest:
    def __init__(self, config_file=False, **kwargs):
        self.studio_app = None
        self.xml_creator = None
        self.bin_path = root_dir
        self.config_path = os.path.join(root_dir, '../config')
        self.path_to_xml = config_file

        self.launch_studio(**kwargs)

        self.create_config_tab_radiobuttons()
        print("Studio launched")

    def launch_studio(self, config_file=False, studio_flag=True, skip_validate_flag=False, rules_flag=True, model3D_flag=False, tensor_flag=False, exec_file="project", nanohub_flag=False, is_movable_flag=False, pytest_flag=False, biwt_flag=False):
        QtCore.QDir.addSearchPath('images', os.path.join(root_dir, 'images'))
        QtCore.QDir.addSearchPath('icon', os.path.join(root_dir, 'icon'))
        self.studio_app = QApplication(sys.argv)
        self.xml_creator = studio.PhysiCellXMLCreator(config_file, studio_flag, skip_validate_flag, rules_flag, model3D_flag, tensor_flag, exec_file, nanohub_flag, is_movable_flag, pytest_flag, biwt_flag)
        # self.studio_app.exec_()
    
    def change_xml(self, filename):
        self.path_to_xml = os.path.join(self.config_path, filename)
        self.xml_creator.current_xml_file = self.path_to_xml
        self.xml_creator.config_file = filename
        if self.xml_creator.studio_flag:
            self.xml_creator.run_tab.config_file = filename
            self.xml_creator.run_tab.config_xml_name.setText(filename)
        self.xml_creator.show_sample_model()
        self.xml_creator.vis_tab.update_output_dir(self.xml_creator.config_tab.folder.text())

    def parse_xml(self):
        tree = ET.parse(self.path_to_xml)
        root = tree.getroot()
        return root

    def check_import(self):
        print("Checking import")
        root = self.parse_xml()
        tags = root.findall('*')
        for tag in tags:
            assert hasattr(self, f"check_{tag.tag}")
            getattr(self, f"check_{tag.tag}")(root)

    def check_domain(self, root):
        domain = root.find('domain')
        for coord in ['x', 'y', 'z']:
            for direction in ['min', 'max']:
                assert domain.find(coord + '_' + direction).text == self.xml_creator.config_tab.__getattribute__(coord + direction).text()
            assert domain.find(f'd{coord}').text == self.xml_creator.config_tab.__getattribute__(f'{coord}del').text()

    def check_overall(self, root):
        overall = root.find('overall')
        assert overall.find('max_time').text == self.xml_creator.config_tab.max_time.text()
        for step in ['diffusion', 'mechanics', 'phenotype']:
            assert overall.find(f"dt_{step}").text == self.xml_creator.config_tab.__getattribute__(f"{step}_dt").text()

    def check_parallel(self, root):
        parallel = root.find('parallel')
        assert parallel.find('omp_num_threads').text == self.xml_creator.config_tab.num_threads.text()

    def check_save(self, root):
        save = root.find('save')
        assert save.find('folder').text == self.xml_creator.config_tab.folder.text()
        full_data = save.find('full_data')
        assert full_data.find('interval').text == self.xml_creator.config_tab.full_interval.text()
        assert (full_data.find('enable').text=='true') == self.xml_creator.config_tab.save_full.isChecked()
        SVG = save.find('SVG')
        if SVG:
            assert SVG.find('interval').text == self.xml_creator.config_tab.svg_interval.text()
            assert (SVG.find('enable').text=='true') == self.xml_creator.config_tab.save_svg.isChecked()
            plot_substrate = SVG.find('plot_substrate')
            if plot_substrate:
                assert (plot_substrate.attrib['enabled']=='true') == self.xml_creator.config_tab.plot_substrate_svg.isChecked()
                assert plot_substrate.find('substrate').text == self.xml_creator.config_tab.svg_substrate_to_plot_dropdown.currentText()
                limits_enabled = plot_substrate.attrib['limits']=='true'
                assert limits_enabled == self.xml_creator.config_tab.plot_substrate_limits.isChecked()
                if limits_enabled:
                    assert plot_substrate.find('min_conc').text == self.xml_creator.config_tab.svg_substrate_min.text()
                    assert plot_substrate.find('max_conc').text == self.xml_creator.config_tab.svg_substrate_max.text()

    def check_options(self, root):
        options = root.find('options')
        assert (options.find('virtual_wall_at_domain_edge').text=='true') == self.xml_creator.config_tab.virtual_walls.isChecked()
        random_seed = options.find('random_seed')
        if random_seed:
            xml_text = random_seed.text
            # check if the seed is an integer
            if xml_text.isdigit():
                assert self.xml_creator.config_tab.random_seed_integer_button.isChecked()
                assert int(xml_text) == int(self.xml_creator.config_tab.random_seed_integer.text())
            else:
                assert self.xml_creator.config_tab.random_seed_random_button.isChecked()

    def check_microenvironment_setup(self, root):
        microenvironment_setup = root.find('microenvironment_setup')
        # get all elements with tag "variable"
        for variable in microenvironment_setup.findall('variable'):
            self.check_microenvironment_variable(variable)
        options = microenvironment_setup.find('options')
        if options:
            self.check_microenvironment_options(options)

    def check_microenvironment_variable(self, variable):
        name = variable.attrib['name']
        physical_parameter_set = variable.find('physical_parameter_set')
        assert physical_parameter_set.find('diffusion_coefficient').text == self.xml_creator.microenv_tab.param_d[name]['diffusion_coef']
        assert physical_parameter_set.find('decay_rate').text == self.xml_creator.microenv_tab.param_d[name]['decay_rate']
        assert variable.find('initial_condition').text == self.xml_creator.microenv_tab.param_d[name]['init_cond']
        Dirichlet_boundary_condition = variable.find('Dirichlet_boundary_condition')
        assert (Dirichlet_boundary_condition.attrib['enabled'].lower()=='true') == self.xml_creator.microenv_tab.param_d[name]['dirichlet_enabled']
        assert Dirichlet_boundary_condition.text == self.xml_creator.microenv_tab.param_d[name]['dirichlet_bc']
        Dirichlet_options = variable.find('Dirichlet_options')
        if Dirichlet_options:
            for boundary_value in Dirichlet_options.findall('boundary_value'):
                id = boundary_value.attrib['ID']
                assert (boundary_value.attrib['enabled'].lower() == 'true') == self.xml_creator.microenv_tab.param_d[name][f'enable_{id}']
                assert boundary_value.text == self.xml_creator.microenv_tab.param_d[name][f'dirichlet_{id}']
            
    def check_microenvironment_options(self, options):
        calculate_gradients = options.find('calculate_gradients')
        if calculate_gradients:
            do_calculate_gradients = calculate_gradients.text.lower() == 'true'
            assert do_calculate_gradients == self.xml_creator.microenv_tab.gradients.isChecked()
            assert do_calculate_gradients == self.xml_creator.microenv_tab.param_d["gradients"]
        track_internalized_substrates_in_each_agent = options.find('track_internalized_substrates_in_each_agent')
        if track_internalized_substrates_in_each_agent:
            do_track_internalized_substrates = track_internalized_substrates_in_each_agent.text.lower() == 'true'
            assert do_track_internalized_substrates == self.xml_creator.microenv_tab.track_in_agents.isChecked()
            assert do_track_internalized_substrates == self.xml_creator.microenv_tab.param_d["track_in_agents"]
        initial_condition = options.find('initial_condition')
        if initial_condition and (initial_condition.attrib['enabled'].lower() == 'true'):
            self.check_microenvironment_initial_condition(initial_condition)
        dirichlet_nodes = options.find('dirichlet_nodes')
        if dirichlet_nodes and (dirichlet_nodes.attrib['enabled'].lower() == 'true'):
            self.check_microenvironment_dirichlet_nodes(dirichlet_nodes)

    def check_microenvironment_initial_condition(self, initial_condition):
        assert self.xml_creator.ics_tab.enable_csv_for_substrate_ics
        assert initial_condition.find('filename').text == self.xml_creator.ics_tab.full_substrate_ic_fname
        # no check for filetype yet?? why am I feeling lazy about this?
    
    def check_microenvironment_dirichlet_nodes(self, dirichlet_nodes):
        # I believe this is not implemented in PhysiCell, so I will not check for it
        raise NotImplementedError
    
    def check_cell_definitions(self, root):
        cell_definitions = root.find('cell_definitions')
        for cell_definition in cell_definitions.findall('cell_definition'):
            self.check_cell_definition(cell_definition)

    def check_cell_definition(self, cell_definition):
        name = cell_definition.attrib['name']
        phenotype = cell_definition.find('phenotype')
        self.check_phenotype(name, phenotype)
        custom_data = cell_definition.find('custom_data')
        if custom_data:
            self.check_custom_data(custom_data, name)
        initial_parameter_distributions = cell_definition.find('initial_parameter_distributions')
        if initial_parameter_distributions:
            self.check_initial_parameter_distributions(initial_parameter_distributions, name)
        
    def check_phenotype(self, name, phenotype):
        tags = ['cycle', 'death', 'volume', 'mechanics', 'motility', 'secretion', 'cell_interactions', 'cell_transformations', 'cell_integrity']
        for tag in tags:
            getattr(self, f"check_{tag}")(name, phenotype)

    def check_cycle(self, name, phenotype):
        cycle_code_dict = {
            0: {'full_name': 'advanced Ki67', 'short_name': 'advancedKi67', 'idx': 2, 'num_phases': 3},
            1: {'full_name': 'basic Ki67', 'short_name': 'Ki67', 'idx': 1, 'num_phases': 2},
            2: {'full_name': 'flow cytometry', 'short_name': 'flowcyto', 'idx': 3, 'num_phases': 3},
            5: {'full_name': 'live cells', 'short_name': 'live', 'idx': 0, 'num_phases': 1},
            6: {'full_name': 'flow cytometry separated', 'short_name': 'flowcytosep', 'idx': 4, 'num_phases': 4},
            7: {'full_name': 'cycling quiescent', 'short_name': 'quiescent', 'idx': 5, 'num_phases': 2}
        } # keys are the PhysiCell codes, the idx is the index of the cycle in the combo box
        cycle = phenotype.find('cycle')
        xml_code = cycle.attrib['code']
        studio_combo_widget_idx = self.xml_creator.celldef_tab.param_d[name]["cycle_choice_idx"]
        studio_code = self.xml_creator.celldef_tab.cycle_combo_idx_code[studio_combo_widget_idx]
        xml_name = cycle.attrib['name']
        studio_name = self.xml_creator.celldef_tab.cycle_combo_idx_name[studio_combo_widget_idx]
        assert xml_code == studio_code
        assert xml_name == studio_name
        is_durations = cycle.find('phase_durations')
        assert (is_durations is not None) == self.xml_creator.celldef_tab.param_d[name]['cycle_duration_flag']
        if is_durations:
            self.check_cycle_durations(name, cycle, cycle_code_dict[int(xml_code)]["short_name"])
        else:
            self.check_cycle_transition_rates(name, cycle, cycle_code_dict[int(xml_code)]["short_name"])

    def check_durations(self, name, model, prefix, cyclic=True):
        phase_durations = model.find('phase_durations')
        durations = phase_durations.findall('duration')
        for duration in durations:
            start_index = int(duration.attrib['index'])
            end_index = (start_index + 1) 
            if cyclic:
                end_index = end_index % len(durations)
            value = duration.text
            cycle_key_base = f'{prefix}_{start_index}{end_index}'
            assert (duration.attrib['fixed_duration'].lower() == 'true') == self.xml_creator.celldef_tab.param_d[name][f'{cycle_key_base}_fixed_duration']
            assert value == self.xml_creator.celldef_tab.param_d[name][f'{cycle_key_base}_duration']

    def check_cycle_durations(self, name, cycle, short_name):
        self.check_durations(name, cycle, f"cycle_{short_name}")
    
    def check_cycle_transition_rates(self, name, cycle, short_name):
        phase_transition_rates = cycle.find('phase_transition_rates')
        rates = phase_transition_rates.findall('rate')
        for rate in rates:
            start_index = int(rate.attrib['start_index'])
            end_index = int(rate.attrib['end_index'])
            value = rate.text
            cycle_key_base = f'cycle_{short_name}_{start_index}{end_index}'
            assert (rate.attrib['fixed_duration'].lower() == 'true') == self.xml_creator.celldef_tab.param_d[name][f'{cycle_key_base}_fixed_trate']
            assert value == self.xml_creator.celldef_tab.param_d[name][f'{cycle_key_base}_trate']
    
    def check_death(self, name, phenotype):
        death = phenotype.find('death')
        for model in death.findall('model'):
            self.check_death_model(name, model)
    
    def check_death_model(self, name, model):
        code = model.attrib['code']
        model_name = model.attrib['name']
        assert model.find('death_rate').text == self.xml_creator.celldef_tab.param_d[name][f'{model_name}_death_rate']
        is_durations = model.find('phase_durations')
        assert (is_durations is not None) == self.xml_creator.celldef_tab.param_d[name][f'{model_name}_duration_flag']
        if is_durations:
            self.check_death_durations(name, model, model_name)
        else:
            self.check_death_transition_rates(name, model, model_name)
        self.check_death_parameters(name, model, model_name)

    def check_death_durations(self, name, model, model_name):
        self.check_durations(name, model, model_name, cyclic=False)

    def check_death_transition_rates(self, name, model, model_name):
        raise NotImplementedError
        
    def check_death_parameters(self, name, model, model_name):
        parameters = model.find('parameters')
        xml_tags = ['unlysed_fluid_change_rate', 'lysed_fluid_change_rate', 'cytoplasmic_biomass_change_rate', 'nuclear_biomass_change_rate', 'calcification_rate', 'relative_rupture_volume']
        key_suffixes = ['unlysed_rate', 'lysed_rate', 'cyto_rate', 'nuclear_rate', 'calcif_rate', 'rel_rupture_volume']
        for xml_tag, key_suffix in zip(xml_tags, key_suffixes):
            assert parameters.find(xml_tag).text == self.xml_creator.celldef_tab.param_d[name][f'{model_name}_{key_suffix}']
        
    def check_volume(self, name, phenotype):
        volume = phenotype.find('volume')
        xml_tags = ['total', 'fluid_fraction', 'nuclear', 'fluid_change_rate', 'cytoplasmic_biomass_change_rate', 'nuclear_biomass_change_rate', 'calcified_fraction', 'calcification_rate', 'relative_rupture_volume']
        keys = ['volume_total', 'volume_fluid_fraction', 'volume_nuclear', 'volume_fluid_change_rate', 'volume_cytoplasmic_rate', 'volume_nuclear_rate', 'volume_calcif_fraction', 'volume_calcif_rate', 'volume_rel_rupture_vol']
        for xml_tag, key in zip(xml_tags, keys):
            assert volume.find(xml_tag).text == self.xml_creator.celldef_tab.param_d[name][key]

    def check_mechanics(self, name, phenotype):
        mechanics = phenotype.find('mechanics')
        xml_tags = ['cell_cell_adhesion_strength', 'cell_cell_repulsion_strength', 'relative_maximum_adhesion_distance', 'attachment_elastic_constant', 'attachment_rate', 'detachment_rate', 'maximum_number_of_attachments']
        keys = ['mechanics_adhesion', 'mechanics_repulsion', 'mechanics_adhesion_distance', 'mechanics_elastic_constant', 'attachment_rate', 'detachment_rate', 'mechanics_max_num_attachments']
        for xml_tag, key in zip(xml_tags, keys):
            assert mechanics.find(xml_tag).text == self.xml_creator.celldef_tab.param_d[name][key]
        cell_adhesion_affinities = mechanics.find('cell_adhesion_affinities')
        if cell_adhesion_affinities:
            self.check_cell_adhesion_affinities(name, cell_adhesion_affinities)
        options = mechanics.find('options')
        if options:
            self.check_mechanics_options(name, options)

    def check_cell_adhesion_affinities(self, name, cell_adhesion_affinities):
        for cell_adhesion_affinity in cell_adhesion_affinities.findall('cell_adhesion_affinity'):
            other_cell_type = cell_adhesion_affinity.attrib['name']
            affinity = cell_adhesion_affinity.text
            assert affinity == self.xml_creator.celldef_tab.param_d[name]['cell_adhesion_affinity'][other_cell_type]

    def check_mechanics_options(self, name, options):
        common_suffixes = ['relative_equilibrium_distance', 'absolute_equilibrium_distance']
        xml_tags = [f"set_{suffix}" for suffix in common_suffixes]
        keys = [f"mechanics_{suffix}" for suffix in common_suffixes]
        enabled_keys = [f"{key}_enabled" for key in keys]
        for xml_tag, key, enabled_key in zip(xml_tags, keys, enabled_keys):
            assert (options.find(xml_tag).attrib['enabled'].lower() == 'true') == self.xml_creator.celldef_tab.param_d[name][enabled_key]
            assert options.find(xml_tag).text == self.xml_creator.celldef_tab.param_d[name][key]
        
    def check_motility(self, name, phenotype):
        motility = phenotype.find('motility')
        xml_tags = ['speed', 'persistence_time', 'migration_bias']
        for tag in xml_tags:
            assert motility.find(tag).text == self.xml_creator.celldef_tab.param_d[name][tag]
        options = motility.find('options')
        if options:
            self.check_motility_options(name, options)
    
    def check_motility_options(self, name, options):
        assert (options.find('enabled').text.lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['motility_enabled']
        assert (options.find('use_2D').text.lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['motility_use_2D']
        chemotaxis = options.find('chemotaxis')
        if chemotaxis:
            self.check_chemotaxis(name, chemotaxis)
        advanced_chemotaxis = options.find('advanced_chemotaxis')
        if advanced_chemotaxis:
            self.check_advanced_chemotaxis(name, advanced_chemotaxis)
    
    def check_chemotaxis(self, name, chemotaxis):
        assert (chemotaxis.find('enabled').text.lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['motility_chemotaxis']
        assert chemotaxis.find('substrate').text == self.xml_creator.celldef_tab.param_d[name]['motility_chemotaxis_substrate']
        assert (chemotaxis.find('direction').text=="1") == self.xml_creator.celldef_tab.param_d[name]['motility_chemotaxis_towards']
    
    def check_advanced_chemotaxis(self, name, advanced_chemotaxis):
        assert (advanced_chemotaxis.find('enabled').text.lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['motility_advanced_chemotaxis']
        assert (advanced_chemotaxis.find('normalize_each_gradient').text.lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['normalize_each_gradient']
        chemotactic_sensitivities = advanced_chemotaxis.find('chemotactic_sensitivities')
        if chemotactic_sensitivities:
            self.check_chemotactic_sensitivities(name, chemotactic_sensitivities)
    
    def check_chemotactic_sensitivities(self, name, chemotactic_sensitivities):
        for chemotactic_sensitivity in chemotactic_sensitivities.findall('chemotactic_sensitivity'):
            substrate = chemotactic_sensitivity.attrib['substrate']
            sensitivity = chemotactic_sensitivity.text
            assert sensitivity == self.xml_creator.celldef_tab.param_d[name]['chemotactic_sensitivity'][substrate]
    
    def check_secretion(self, name, phenotype):
        secretion = phenotype.find('secretion')
        for substrate in secretion.findall('substrate'):
            self.check_secretion_substrate(name, substrate)
    
    def check_secretion_substrate(self, name, substrate):
        substrate_name = substrate.attrib['name']
        tags = ['secretion_rate', 'secretion_target', 'uptake_rate', 'net_export_rate']
        for tag in tags:
            assert substrate.find(tag).text == self.xml_creator.celldef_tab.param_d[name]['secretion'][substrate_name][tag]
    
    def check_cell_interactions(self, name, phenotype):
        cell_interactions = phenotype.find('cell_interactions')
        tags = ['apoptotic_phagocytosis_rate', 'necrotic_phagocytosis_rate', 'other_dead_phagocytosis_rate', 'attack_damage_rate', 'attack_duration']
        for tag in tags:
            assert cell_interactions.find(tag).text == self.xml_creator.celldef_tab.param_d[name][tag]
        for block in ['live_phagocytosis_rates', 'attack_rates', 'fusion_rates']:
            self.check_cell_interactions_rates_block(name, cell_interactions, block)
    
    def check_cell_interactions_rates_block(self, name, cell_interactions, tag):
        rates_block = cell_interactions.find(tag)
        if not rates_block:
            return
        key = tag.rstrip('s')
        for rate in rates_block.findall('*'):
            other_cell_type = rate.attrib['name']
            rate_value = rate.text
            assert rate_value == self.xml_creator.celldef_tab.param_d[name][key][other_cell_type]

    def check_cell_transformations(self, name, phenotype):
        cell_transformations = phenotype.find('cell_transformations')
        self.check_cell_interactions_rates_block(name, cell_transformations, 'transformation_rates')
    
    def check_cell_integrity(self, name, phenotype):
        cell_integrity = phenotype.find('cell_integrity')
        tags = ['damage_rate', 'damage_repair_rate']
        for tag in tags:
            assert cell_integrity.find(tag).text == self.xml_creator.celldef_tab.param_d[name][tag]
    
    def check_custom_data(self, custom_data, name):
        for tag in custom_data.findall('*'):
            assert tag.text == self.xml_creator.celldef_tab.param_d[name]['custom_data'][tag.tag][0]
            assert (tag.attrib['conserved'].lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['custom_data'][tag.tag][1]
    
    def check_initial_parameter_distributions(self, initial_parameter_distributions, name):
        assert (initial_parameter_distributions.attrib['enabled'].lower() == 'false') == self.xml_creator.celldef_tab.param_d[name]['par_dists_disabled']
        for distribution in initial_parameter_distributions.findall('distribution'):
            self.check_initial_parameter_distribution(name, distribution)

    def check_initial_parameter_distribution(self, name, distribution):
        dist_dict = {"uniform": "Uniform", "loguniform": "Log Uniform", "normal": "Normal", "lognormal": "Log Normal", "log10normal": "Log10 Normal"}
        behavior = distribution.find('behavior').text
        assert (distribution.attrib['enabled'].lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['par_dists'][behavior]['enabled']
        assert dist_dict[distribution.attrib['type'].lower()] == self.xml_creator.celldef_tab.param_d[name]['par_dists'][behavior]['distribution']
        assert (distribution.attrib['check_base'].lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['par_dists'][behavior]['enforce_base']
        for tag in distribution.findall('*'):
            if tag.tag == 'behavior':
                continue
            assert tag.text == self.xml_creator.celldef_tab.param_d[name]['par_dists'][behavior]['parameters'][tag.tag]

    def check_initial_conditions(self, root):
        initial_conditions = root.find('initial_conditions')
        cell_positions = initial_conditions.find('cell_positions')
        assert (cell_positions.attrib['enabled'].lower() == 'true') == self.xml_creator.config_tab.cells_csv.isChecked()
        assert cell_positions.find('folder').text == self.xml_creator.config_tab.csv_folder.text()
        assert cell_positions.find('filename').text == self.xml_creator.config_tab.csv_file.text()
    
    def check_cell_rules(self, root):
        cell_rules = root.find('cell_rules')
        rulesets = cell_rules.find('rulesets')
        ruleset_s = rulesets.findall('ruleset')
        assert len(ruleset_s) == 1
        ruleset = ruleset_s[0]
        assert (ruleset.attrib['enabled'].lower() == 'true') == self.xml_creator.rules_tab.rules_enabled.isChecked()
        assert (ruleset.attrib['enabled'].lower() == 'true') == self.xml_creator.rules_tab.rules_enabled_attr
        assert ruleset.find('folder').text == self.xml_creator.rules_tab.rules_folder.text()
        assert ruleset.find('filename').text == self.xml_creator.rules_tab.rules_file.text()
    
    def check_user_parameters(self, root):
        user_parameters = root.find('user_parameters')
        for user_parameter in user_parameters.findall('*'):
            self.check_user_parameter(user_parameter)

    def check_user_parameter(self, user_parameter):
        type = user_parameter.attrib['type']
        value = user_parameter.text
        for irow in range(self.xml_creator.user_params_tab.max_rows):
            if user_parameter.tag == self.xml_creator.user_params_tab.utable.cellWidget(irow, self.xml_creator.user_params_tab.var_icol_name).text():
                assert type == self.xml_creator.user_params_tab.utable.cellWidget(irow, self.xml_creator.user_params_tab.var_icol_type).currentText()
                assert value == self.xml_creator.user_params_tab.utable.cellWidget(irow, self.xml_creator.user_params_tab.var_icol_value).text()
                return
    
    ########### test changes in gui and check xml ############

    def save_xml(self):
        self.xml_creator.celldef_tab.check_valid_cell_defs()
        filepath = os.path.join(self.config_path, '__test__.xml')

        self.xml_creator.celldef_tab.config_path = self.xml_creator.current_xml_file
        self.xml_creator.config_tab.fill_xml()
        self.xml_creator.microenv_tab.fill_xml()
        self.xml_creator.celldef_tab.fill_xml()
        self.xml_creator.user_params_tab.fill_xml()
        if self.xml_creator.rules_flag:
            self.xml_creator.rules_tab.fill_xml()

        self.xml_creator.tree.write(filepath)
        pretty_print(filepath, filepath)

    def compare_xml(self, xml1, xml2, values_to_check):
        tree1 = ET.parse(xml1)
        tree2 = ET.parse(xml2)
        root1 = tree1.getroot()
        root2 = tree2.getroot()
        return self.compare_elements(root1, root2, values_to_check)
    
    def compare_elements(self, elem1, elem2, values_to_check, attributes_to_check={}):
        print(f"Comparing {elem1.tag} and {elem2.tag}")
        if elem1.tag != elem2.tag:
            return False
        t1 = elem1.text
        t2 = elem2.text
        t1_is_nothing = t1 is None or t1.strip() == ''
        t2_is_nothing = t2 is None or t2.strip() == ''
        if (t1_is_nothing != t2_is_nothing) or (t1 != t2):
            print(f"Content does not match: {t1} != {t2}")
            return False
        if elem1.attrib != elem2.attrib:
            # make sure the dicts have the same keys
            if set(elem1.attrib.keys()) != set(elem2.attrib.keys()):
                print(f"Attribute keys do not match: {elem1.attrib.keys()} != {elem2.attrib.keys()}")
                return False
            print(f"Attributes do not match: {elem1.attrib} != {elem2.attrib}")
            print(f"Attributes to ignore: {attributes_to_check}")
            for a in elem1.attrib.keys():
                if (a in attributes_to_check.keys()):
                    if (elem2.attrib[a] != attributes_to_check[a]):
                        print("Did not change attribute correctly")
                        return False
                elif elem1.attrib[a] != elem2.attrib[a]:
                    print("Attribute values do not match")
                    return False 
        if 'enabled' in elem1.attrib.keys() and elem1.attrib['enabled'].lower() == 'false':
            return True # don't bother comparing things that are in disabled blocks
        if len(elem1) != len(elem2):
            print(f"Lengths do not match: {len(elem1)} != {len(elem2)}")
            return False
        paths_to_check = list(values_to_check.keys())
        for child1, child2 in zip(elem1, elem2):
            print(f"{values_to_check = }")
            if child1.tag in values_to_check.keys():
                index = paths_to_check.index(child1.tag)
                print(f"Checking {child1.tag} for {values_to_check[child1.tag]}")
                print(f"Value in xml2: {child2.text}, expected value: {values_to_check[child1.tag]}")
                assert child2.text == values_to_check[child1.tag]
                continue
            temp = [k.startswith(f'{child1.tag}/@') for k in paths_to_check]
            _attributes_to_check = {}
            if any(temp):
                print(f"Checking {child1.tag} for {values_to_check}")
                index = temp.index(True)
                attrib_name = paths_to_check[index].split('/@')[1]
                print(f"Value in xml2: {child2.text}, expected value: {values_to_check}")
                _attributes_to_check[attrib_name] = values_to_check[paths_to_check[index]]
            _values_to_check = {'//'.join(k.split('//')[1:]): v for k, v in values_to_check.items() if k.startswith(child1.tag)}
            print(f"{_values_to_check = }")
            if not self.compare_elements(child1, child2, _values_to_check, attributes_to_check=_attributes_to_check):
                return False
        return True

    ########### perturb in gui and check xml ############
    config_tab_text_fields = {'xmin': 'domain//x_min', 'xmax': 'domain//x_max', 'ymin': 'domain//y_min', 'ymax': 'domain//y_max', 'zmin': 'domain//z_min', 'zmax': 'domain//z_max', 'xdel': 'domain//dx', 'ydel': 'domain//dy', 'zdel': 'domain//dz',
                              'max_time': 'overall//max_time', 'diffusion_dt': 'overall//dt_diffusion', 'mechanics_dt': 'overall//dt_mechanics', 'phenotype_dt': 'overall//dt_phenotype',
                              'num_threads': 'parallel//omp_num_threads', 'folder': 'save//folder', 'full_interval': 'save//full_data//interval', 'svg_interval': 'save//SVG//interval',
                              'svg_substrate_min': 'save//SVG//plot_substrate//min_conc', 'svg_substrate_max': 'save//SVG//plot_substrate//max_conc', 'random_seed_integer': 'options//random_seed',
                              'csv_folder': 'initial_conditions//cell_positions//folder', 'csv_file': 'initial_conditions//cell_positions//filename'}
    config_tab_text_fields = {k: {'xml_path': v, 'new_value': "12321"} for k, v in config_tab_text_fields.items()}
    config_tab_text_fields["zmin"]["new_value"] = "0.12321"
    config_tab_text_fields["zmax"]["new_value"] = "0.12321"
    config_tab_text_fields["full_interval"]["extra_values_to_check"] = {"save//SVG//interval": "12321"}
    config_tab_text_fields["svg_interval"]["extra_values_to_check"] = {"save//full_data//interval": "12321"}

    config_tab_checkboxes = {'save_full': 'save//full_data//enable', 'save_svg': 'save//SVG//enable', 
                             'plot_substrate_svg': 'save//SVG//plot_substrate/@enabled', 
                             'plot_substrate_limits': 'save//SVG//plot_substrate/@limits', 
                             'virtual_walls': 'options//virtual_wall_at_domain_edge', 'cells_csv': 'initial_conditions//cell_positions/@enabled'}
    config_tab_checkboxes = {k: {'xml_path': v, 'new_value': lambda original_value: str(not original_value).lower()} for k, v in config_tab_checkboxes.items()}
    
    def create_config_tab_radiobuttons(self):
        self.config_tab_radiobuttons = {'random_seed_integer_button': 'options//random_seed', 'random_seed_random_button': 'options//random_seed'}
        self.config_tab_radiobuttons = {k: {'xml_path': v} for k, v in self.config_tab_radiobuttons.items()}
        self.config_tab_radiobuttons["random_seed_integer_button"]["new_value"] = lambda original_value: "system_clock" if original_value else "12321"
        self.config_tab_radiobuttons["random_seed_random_button"]["new_value"] = lambda original_value: "12321" if original_value else "system_clock"
        self.config_tab_radiobuttons["random_seed_integer_button"]["toggle_fn"] = self.toggle_random_seed_gp
        self.config_tab_radiobuttons["random_seed_random_button"]["toggle_fn"] = self.toggle_random_seed_gp

    def toggle_random_seed_gp(self):
        current_id = self.xml_creator.config_tab.random_seed_gp.checkedId()
        self.xml_creator.config_tab.random_seed_gp.button((current_id+1) % 2).click()

    config_tab_combo_boxes = ['svg_substrate_to_plot_dropdown']
    
    def __del__(self):
        self.studio_app.quit()

if __name__ == "__main__":
    studio_test = StudioTest(config_file="./config/PhysiCell_settings.xml")
    config_tab = studio_test.xml_creator.config_tab
    studio_test.check_import()
    # search through this file for all 
    config_tab = studio_test.xml_creator.config_tab
    for field, d in studio_test.config_tab_text_fields.items():
        print(f"Checking {field}")
        original_value = config_tab.__getattribute__(field).text()
        config_tab.__getattribute__(field).setText(d['new_value'])
        studio_test.save_xml()
        values_to_check = {d['xml_path']: d['new_value']}
        if 'extra_values_to_check' in d.keys():
            values_to_check.update(d['extra_values_to_check'])
        assert studio_test.compare_xml(studio_test.path_to_xml, "./config/__test__.xml", values_to_check)
        config_tab.__getattribute__(field).setText(original_value)

    for field, d in studio_test.config_tab_checkboxes.items():
        print(f"Checking {field}")
        original_value = config_tab.__getattribute__(field).isChecked()
        config_tab.__getattribute__(field).setChecked(not original_value)
        studio_test.save_xml()
        values_to_check = {d['xml_path']: d['new_value'](original_value)}
        assert studio_test.compare_xml(studio_test.path_to_xml, "./config/__test__.xml", values_to_check)
        config_tab.__getattribute__(field).setChecked(original_value)

    for field, d in studio_test.config_tab_radiobuttons.items():
        print(f"Checking {field}")
        original_value = config_tab.__getattribute__(field).isChecked()
        d["toggle_fn"]()
        studio_test.save_xml()
        values_to_check = {d['xml_path']: d['new_value'](original_value)}
        assert studio_test.compare_xml(studio_test.path_to_xml, "./config/__test__.xml", values_to_check)
        d["toggle_fn"]()

