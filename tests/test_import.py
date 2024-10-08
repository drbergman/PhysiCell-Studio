import pytest

import os
import sys
import contextlib

import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html

# Add the directory containing the module to sys.path
root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin'))
# root_dir = '/Users/danielbergman/PhysiCell-Studio_git/bin'
sys.path.append(root_dir)

# Now you can import the module
import studio
from pretty_print_xml import pretty_print


from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication

@contextlib.contextmanager
def suppress_output():
    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            yield
            
class XMLPathValue:
    def __init__(self, xml_path, value):
        if isinstance(xml_path, str):
            self.str = xml_path
            self.list = self.tokenize_path(xml_path)
        else:
            self.str = '//'.join([':'.join(e) for e in xml_path])
            self.list = xml_path
        self.value = value

    def tokenize_path(self, path):
        # 'cell_definitions//cell_definition:name:default' -> [ ['cell_definitions'], ['cell_definition', 'name', 'default'] ]
        l = path.split('//')
        ret_val = [self.tokenize_element(e) for e in l]
        return ret_val
    
    def tokenize_element(self, element):
        # 'cell_definition:name:default' -> ['cell_definition', 'name', 'default']
        return element.split(':')

    def is_attrib(self):
        return len(self.list[-1]) != 1

    def is_next_level(self, tag, attrib):
        # entering this, self.list = [[<current xml level>], [<next xml level>], ...]
        # we already know the current xml level worked (that's why we just checked it)
        # now see if we want to keep checking this at the next level down
        if len(self.list) == 1:
            # we are done recursing through this path, don't use it deeper
            return False
        if self.list[1][0] != tag:
            # this xml_path does not start with the correct tag, so it is not part of the next level
            return False
        if len(self.list[1]) == 1:
            # this xml_path does not have an attrib, so it is part of the next level
            return True
        if self.list[1][1] not in attrib.keys():
            # this xml_path does not have the correct attrib, so it is not part of the next level
            return False
        if len(self.list[1]) == 2:
            # this xml_path has the correct attrib, but not value, so we will check the value at the next level
            return True
        if attrib[self.list[1][1]] == self.list[1][2]:
            # this xml_path has the correct attrib value, so it is part of the next level
            return True
        return False

    def length(self):
        return len(self.list)

    def next_attrib_value(self):
        if len(self.list[0]) < 3:
            return None
        return self.list[0][2]

class StudioTest:
    def __init__(self, config_file=False, **kwargs):
        self.studio_app = None
        self.xml_creator = None
        self.bin_path = root_dir
        self.config_path = os.path.join(root_dir, '../config')
        self.path_to_xml = config_file

        self.launch_studio(**kwargs)

        self.create_config_tab_widget_dicts()
        self.create_celldef_tab_widget_dicts()
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
        if SVG is not None:
            assert SVG.find('interval').text == self.xml_creator.config_tab.svg_interval.text()
            assert (SVG.find('enable').text=='true') == self.xml_creator.config_tab.save_svg.isChecked()
            plot_substrate = SVG.find('plot_substrate')
            if plot_substrate is not None:
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
        if random_seed is not None:
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
        if options is not None:
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
        if Dirichlet_options is not None:
            for boundary_value in Dirichlet_options.findall('boundary_value'):
                id = boundary_value.attrib['ID']
                assert (boundary_value.attrib['enabled'].lower() == 'true') == self.xml_creator.microenv_tab.param_d[name][f'enable_{id}']
                assert boundary_value.text == self.xml_creator.microenv_tab.param_d[name][f'dirichlet_{id}']
            
    def check_microenvironment_options(self, options):
        calculate_gradients = options.find('calculate_gradients')
        if calculate_gradients is not None:
            do_calculate_gradients = calculate_gradients.text.lower() == 'true'
            assert do_calculate_gradients == self.xml_creator.microenv_tab.gradients.isChecked()
            assert do_calculate_gradients == self.xml_creator.microenv_tab.param_d["gradients"]
        track_internalized_substrates_in_each_agent = options.find('track_internalized_substrates_in_each_agent')
        if track_internalized_substrates_in_each_agent is not None:
            do_track_internalized_substrates = track_internalized_substrates_in_each_agent.text.lower() == 'true'
            assert do_track_internalized_substrates == self.xml_creator.microenv_tab.track_in_agents.isChecked()
            assert do_track_internalized_substrates == self.xml_creator.microenv_tab.param_d["track_in_agents"]
        initial_condition = options.find('initial_condition')
        if (initial_condition is not None) and (initial_condition.attrib['enabled'].lower() == 'true'):
            self.check_microenvironment_initial_condition(initial_condition)
        dirichlet_nodes = options.find('dirichlet_nodes')
        if (dirichlet_nodes is not None) and (dirichlet_nodes.attrib['enabled'].lower() == 'true'):
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
        if custom_data is not None:
            self.check_custom_data(custom_data, name)
        initial_parameter_distributions = cell_definition.find('initial_parameter_distributions')
        if initial_parameter_distributions is not None:
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
        if is_durations is not None:
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

    def check_transition_rates(self, name, model, prefix):
        phase_transition_rates = model.find('phase_transition_rates')
        rates = phase_transition_rates.findall('rate')
        for rate in rates:
            start_index = int(rate.attrib['start_index'])
            end_index = int(rate.attrib['end_index'])
            value = rate.text
            key_base = f'{prefix}_{start_index}{end_index}'
            assert (rate.attrib['fixed_duration'].lower() == 'true') == self.xml_creator.celldef_tab.param_d[name][f'{key_base}_fixed_trate']
            assert float(value) == float(self.xml_creator.celldef_tab.param_d[name][f'{key_base}_trate'])

    def check_cycle_durations(self, name, cycle, short_name):
        self.check_durations(name, cycle, f"cycle_{short_name}")
    
    def check_cycle_transition_rates(self, name, cycle, short_name):
        self.check_transition_rates(name, cycle, f"cycle_{short_name}")
    
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
        if is_durations is not None:
            self.check_death_durations(name, model, model_name)
        else:
            self.check_death_transition_rates(name, model, model_name)
        self.check_death_parameters(name, model, model_name)

    def check_death_durations(self, name, model, model_name):
        self.check_durations(name, model, model_name, cyclic=False)

    def check_death_transition_rates(self, name, model, model_name):
        self.check_transition_rates(name, model, model_name)
        
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
        if cell_adhesion_affinities is not None:
            self.check_cell_adhesion_affinities(name, cell_adhesion_affinities)
        options = mechanics.find('options')
        if options is not None:
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
        if options is not None:
            self.check_motility_options(name, options)
    
    def check_motility_options(self, name, options):
        assert (options.find('enabled').text.lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['motility_enabled']
        assert (options.find('use_2D').text.lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['motility_use_2D']
        chemotaxis = options.find('chemotaxis')
        if chemotaxis is not None:
            self.check_chemotaxis(name, chemotaxis)
        advanced_chemotaxis = options.find('advanced_chemotaxis')
        if advanced_chemotaxis is not None:
            self.check_advanced_chemotaxis(name, advanced_chemotaxis)
    
    def check_chemotaxis(self, name, chemotaxis):
        assert (chemotaxis.find('enabled').text.lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['motility_chemotaxis']
        assert chemotaxis.find('substrate').text == self.xml_creator.celldef_tab.param_d[name]['motility_chemotaxis_substrate']
        assert (chemotaxis.find('direction').text=="1") == self.xml_creator.celldef_tab.param_d[name]['motility_chemotaxis_towards']
    
    def check_advanced_chemotaxis(self, name, advanced_chemotaxis):
        assert (advanced_chemotaxis.find('enabled').text.lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['motility_advanced_chemotaxis']
        assert (advanced_chemotaxis.find('normalize_each_gradient').text.lower() == 'true') == self.xml_creator.celldef_tab.param_d[name]['normalize_each_gradient']
        chemotactic_sensitivities = advanced_chemotaxis.find('chemotactic_sensitivities')
        if chemotactic_sensitivities is not None:
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
        if rates_block is None:
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
        with suppress_output():
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

    def compare_xml(self, xml1, xml2, content_to_compare, attrib_to_compare, terminal_paths, nonexistant_paths=[], rename_dict={}, show_log=False):
        tree1 = ET.parse(xml1)
        tree2 = ET.parse(xml2)
        root1 = tree1.getroot()
        root2 = tree2.getroot()
        ret_val = self.compare_elements(root1, root2, content_to_compare, attrib_to_compare, terminal_paths, nonexistant_paths, rename_dict, show_log)
        if len(nonexistant_paths) > 0:
            print(f"\tSuccessfully did not find the following paths:")
            for nep in nonexistant_paths:
                print(f"\t\t{nep.str}")
        return ret_val

    def compare_elements(self, elem1, elem2, content_to_compare, attrib_to_compare, terminal_paths, nonexistant_paths, rename_dict, show_log):
        try:
            assert all([nep.length()>1 for nep in nonexistant_paths]), f"We reach a path that should be nonexistant at {elem2.tag}"
        except AssertionError as e:
            raise e
        if show_log:
            # print(f"{elem2.tag = }, {elem2.attrib = }")
            if len(terminal_paths) > 0:
                print(f"{[a.list for a in attrib_to_compare] =}")
                print(f"Terminal paths: {[t.list for t in terminal_paths]}")
        self.compare_content(elem1, elem2, content_to_compare, show_log)
        self.compare_attrib(elem1, elem2, attrib_to_compare, rename_dict, show_log)
        if 'enabled' in elem2.attrib.keys() and elem2.attrib['enabled'].lower() == 'false':
            if show_log:
                print(f"Element {elem2.tag} is disabled")
            return True # if the new element is disabled, once we have checked content and attribs, we can return True
        if len(terminal_paths)==1 and terminal_paths[0].length()==1:
            # this path is terminal if and only if there is a single terminal path and it has one element left t
            print(f"\tTerminating recursion at {':'.join(terminal_paths[0].list[0])}")
            return True
        if show_log and elem1.tag == 'cell_definitions':
            print(f"\tElem1 cell definitions = {[c.attrib['name'] for c in elem1]}")
            print(f"\tElem2 cell definitions = {[c.attrib['name'] for c in elem2]}")
            c1, c2 = self.zip_children(elem1, elem2)
            print(f"\t{[c.attrib['name'] for c in c1] = }")
            print(f"\t{[c.attrib['name'] for c in c2] = }")
        for child1, child2 in self.zip_children(elem1, elem2):
            _content_to_compare = [XMLPathValue(p.list[1:], p.value) for p in content_to_compare if p.is_next_level(child2.tag, child2.attrib)]
            _attrib_to_compare = [XMLPathValue(p.list[1:], p.value) for p in attrib_to_compare if p.is_next_level(child2.tag, child2.attrib)]
            _terminal_paths = [XMLPathValue(p.list[1:], None) for p in terminal_paths if p.is_next_level(child2.tag, child2.attrib)]
            _nonexistant_paths = [XMLPathValue(p.list[1:], None) for p in nonexistant_paths if p.is_next_level(child2.tag, child2.attrib)]
            assert self.compare_elements(child1, child2, _content_to_compare, _attrib_to_compare, _terminal_paths, _nonexistant_paths, rename_dict, show_log)

        return True

    def zip_children(self, elem1, elem2):
        if len(elem2) == 0: # if no children in new xml, then return an empty iterator
            return zip([], [])
        if 'enabled' in elem1.attrib.keys() and elem1.attrib['enabled'].lower() == 'false':
            print(f"\t{elem1.tag} is disabled in the old xml. Recursing through new xml to check if values and attribs match those explicitly set.")
            return zip(elem2, elem2)
        if 'name' not in elem2[0].attrib.keys():
            return zip(elem1, elem2)
        elem2_children = elem2.findall('*')
        elem1_children = []
        for child in elem2_children:
            if child.attrib == {}: # no attribs to compare
                elem1_child = elem1.find(child.tag)
            elif 'name' in child.attrib.keys():
                elem1_child = elem1.find(f".//{child.tag}[@name='{child.attrib['name']}']")
            elif 'ID' in child.attrib.keys():
                elem1_child = elem1.find(f".//{child.tag}[@ID='{child.attrib['ID']}']")
            elif 'index' in child.attrib.keys():
                elem1_child = elem1.find(f".//{child.tag}[@index='{child.attrib['index']}']")
            else:
                elem1_child = elem1.find(child.tag)
            if elem1_child is not None:
                elem1_children.append(elem1_child)
            else:
                # if no element found in original xml, then "compare" the new one with itself. it should terminate the check
                elem1_children.append(child)
        return zip(elem1_children, elem2_children)

    def compare_attrib(self, elem1, elem2, attrib_to_compare: list[XMLPathValue], rename_dict, show_log):
        # if show_log:
        #     print(f"Comparing attrib: {xml_path_value.list}")
        for xml_path_value in attrib_to_compare:
            if show_log:
                print(f"\t{xml_path_value.list = }")
                print(f"\t{xml_path_value.next_attrib_value() = }")
            if xml_path_value.next_attrib_value() in rename_dict.keys():
                # on new cell, the attrib of name in cell_definition will not match the "new value" because it is for the next cell def
                continue
            if xml_path_value.length()==1:
                if show_log:
                    print(f"\t{elem1.attrib = }, \n\t{elem2.attrib = }, \n\t{[a.list for a in attrib_to_compare] = }")
                print(f"\tAttribute assertion: {elem2.attrib[xml_path_value.list[0][1]]} == {xml_path_value.value}")
                try:
                    assert elem2.attrib[xml_path_value.list[0][1]] == xml_path_value.value, f"Attrib {elem2.tag}[@{xml_path_value.list[0][1]}] does not match expected value: {elem2.attrib[xml_path_value.list[0][1]]} != {xml_path_value.value}"
                except AssertionError as e:
                    raise e
                return
        # otherwise, compare the two elem attribs
        differs = False
        # make sure that all elem2 attributes are in elem1
        differs = any([k not in elem1.attrib.keys() for k in elem2.attrib.keys()])
        if not differs:
            for key in elem1.attrib.keys():
                if key not in elem2.attrib.keys():
                    differs = True
                    break
                if elem1.attrib[key] != elem2.attrib[key]:
                    # make sure attrib is not renamed
                    if elem1.attrib[key] in rename_dict.keys() and rename_dict[elem1.attrib[key]] == elem2.attrib[key]:
                        continue
                    differs = True
                    break
        try:
            assert not differs, f"Attribs do not match for {elem1.tag}: {elem1.attrib} != {elem2.attrib}"
        except AssertionError as e:
            print(f"Attributes to compare: {attrib_to_compare}")
            raise e
        return

    def compare_content(self, elem1, elem2, content_to_compare: list[XMLPathValue], show_log):
        for xml_path_value in content_to_compare:
            if xml_path_value.length()==1:
                print(f"\tContent assertion: {elem2.text} == {xml_path_value.value}")
                assert elem2.text == xml_path_value.value
                return
        # otherwise, compare the two elem contents       
        differs = False
        t1 = elem1.text
        t2 = elem2.text
        t1_is_nothing = t1 is None or t1.strip() == ''
        t2_is_nothing = t2 is None or t2.strip() == ''
        if (t1_is_nothing != t2_is_nothing) or (t1 != t2):
            # check if t1 and t2 are numbers
            try:
                if float(t1) != float(t2):
                    differs = True
            except ValueError:
                differs = True
        try:
            assert not differs, f"Content does not match: {t1} != {t2}"
        except AssertionError as e:
            raise e
        return 
    
    def tokenize_path(self, path):
        return path.split('//')
    
    def tokenize_element(self, element):
        return element.split(':')
    
    def test_tab(self, tab_name, tab_path: list[str]):
        print(f"\n{'-'*10}TESTING {tab_name.upper()}{'-'*10}\n")
        tab = self.xml_creator
        for p in tab_path:
            tab = tab.__getattribute__(p)
        for suffix in ['text_fields', 'checkboxes', 'radiobuttons', 'comboboxes', 'trees', 'buttons']:
            if hasattr(self, f"{tab_name}_{suffix}"):
                testing_fn = getattr(self, f"test_{suffix}")
                print(f"\n\t-- Testing {tab_name}_{suffix}:\n")
                items_to_test = getattr(self, f"{tab_name}_{suffix}")
                testing_fn(items_to_test, tab)
            else:
                print(f"No {tab_name}_{suffix} to test")

    def test_text_fields(self, items_to_test, tab):
        for field, d in items_to_test.items():
            print(f"Checking {field}")
            original_value = tab.__getattribute__(field).text()
            if 'pre_set_fn' in d.keys():
                d['pre_set_fn']()
            tab.__getattribute__(field).setText(d['new_value'])
            self.save_xml()
            xml_path_value = XMLPathValue('PhysiCell_settings//' + d['xml_path'], d['new_value'])
            if xml_path_value.is_attrib():
                attrib_to_compare = [xml_path_value]
                content_to_compare = []
            else:
                content_to_compare = [xml_path_value]
                attrib_to_compare = []
            if 'extra_content_to_compare' in d.keys():
                content_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_content_to_compare'].items()]
            if 'extra_attrib_to_compare' in d.keys():
                attrib_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_attrib_to_compare'].items()]
            terminal_paths = []
            assert self.compare_xml(self.path_to_xml, "./config/__test__.xml", content_to_compare, attrib_to_compare, terminal_paths, nonexistant_paths=[], rename_dict={}, show_log=False)
            if 'post_check_fn' in d.keys():
                d['post_check_fn']()
            tab.__getattribute__(field).setText(original_value)
            
    def test_checkboxes(self, items_to_test, tab):
        for field, d in items_to_test.items():
            print(f"Checking {field}")
            original_value = tab.__getattribute__(field).isChecked()
            tab.__getattribute__(field).setChecked(not original_value)
            self.save_xml()
            xml_path_value = XMLPathValue('PhysiCell_settings//' + d['xml_path'], d['new_value'](original_value))
            if xml_path_value.is_attrib():
                attrib_to_compare = [xml_path_value]
                content_to_compare = []
            else:
                content_to_compare = [xml_path_value]
                attrib_to_compare = []
            if 'extra_content_to_compare' in d.keys():
                content_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_content_to_compare'].items()]
            if 'extra_attrib_to_compare' in d.keys():
                attrib_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_attrib_to_compare'].items()]
            assert self.compare_xml(self.path_to_xml, "./config/__test__.xml", content_to_compare, attrib_to_compare, [])
            tab.__getattribute__(field).setChecked(original_value)

    def test_radiobuttons(self, items_to_test, tab):
        for field, d in items_to_test.items():
            print(f"Checking {field}")
            original_value = tab.__getattribute__(field).isChecked()
            d["toggle_fn"]()
            self.save_xml()
            xml_path_value = XMLPathValue('PhysiCell_settings//' + d['xml_path'], d['new_value'](original_value))
            if xml_path_value.is_attrib():
                attrib_to_compare = [xml_path_value]
                content_to_compare = []
            else:
                content_to_compare = [xml_path_value]
                attrib_to_compare = []
            if 'extra_content_to_compare' in d.keys():
                content_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_content_to_compare'].items()]
            if 'extra_attrib_to_compare' in d.keys():
                attrib_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_attrib_to_compare'].items()]
            assert self.compare_xml(self.path_to_xml, "./config/__test__.xml", content_to_compare, attrib_to_compare, [])
            d["toggle_fn"]()

    def test_comboboxes(self, items_to_test, tab):
        for field, d in items_to_test.items():
            print(f"Checking {field}")
            original_value = tab.__getattribute__(field).currentText()
            original_index = tab.__getattribute__(field).currentIndex()
            if 'pre_loop_fn' in d.keys():
                d['pre_loop_fn']()
            for i in range(tab.__getattribute__(field).count()):
                tab.__getattribute__(field).setCurrentIndex((original_index + i) % tab.__getattribute__(field).count())
                new_value = tab.__getattribute__(field).currentText()
                self.save_xml()
                xml_path_value = XMLPathValue('PhysiCell_settings//' + d['xml_path'], new_value)
                if xml_path_value.is_attrib():
                    attrib_to_compare = [xml_path_value]
                    content_to_compare = []
                else:
                    content_to_compare = [xml_path_value]
                    attrib_to_compare = []
                if 'extra_content_to_compare' in d.keys():
                    content_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_content_to_compare'].items()]
                if 'extra_attrib_to_compare' in d.keys():
                    attrib_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_attrib_to_compare'].items()]
                assert self.compare_xml(self.path_to_xml, "./config/__test__.xml", content_to_compare, attrib_to_compare, [], show_log=False)
            tab.__getattribute__(field).setCurrentText(original_value)
            tab.__getattribute__(field).setCurrentIndex(original_index)
            if 'post_loop_fn' in d.keys():
                d['post_loop_fn']()

    def test_trees(self, items_to_test, tab):
        for field, d in items_to_test.items():
            print(f"Checking {field}")
            tree = tab.__getattribute__(field)
            for i in range(tab.__getattribute__(field).topLevelItemCount()):
                original_value = tree.topLevelItem(i).text(0)
                new_value = f"{original_value}__test__"
                with suppress_output():
                    tree.topLevelItem(i).setText(0, new_value)
                self.save_xml()
                xml_path_value = XMLPathValue('PhysiCell_settings//' + d['xml_path'], new_value)
                if xml_path_value.is_attrib():
                    attrib_to_compare = [xml_path_value]
                    content_to_compare = []
                else:
                    content_to_compare = [xml_path_value]
                    attrib_to_compare = []
                if 'extra_content_to_compare' in d.keys():
                    content_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_content_to_compare'].items()]
                if 'extra_attrib_to_compare' in d.keys():
                    attrib_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_attrib_to_compare'].items()]
                rename_dict = {original_value: new_value}
                terminal_paths = []
                assert self.compare_xml(self.path_to_xml, "./config/__test__.xml", content_to_compare, attrib_to_compare, terminal_paths, nonexistant_paths=[], rename_dict=rename_dict, show_log=False)
                with suppress_output():
                    tree.topLevelItem(i).setText(0, original_value)

    def test_buttons(self, items_to_test, tab):
        for field, d in items_to_test.items():
            print(f"Checking {field}")
            with suppress_output():
                tab.__getattribute__(field).click()
            if 'post_push_fn' in d.keys():
                with suppress_output():
                    d['post_push_fn']()
            self.save_xml()
            content_to_compare = []
            attrib_to_compare = []
            if 'xml_path' in d.keys():
                xml_path_value = XMLPathValue('PhysiCell_settings//' + d['xml_path'], d['new_value'])
                if xml_path_value.is_attrib():
                    attrib_to_compare = [xml_path_value]
                else:
                    content_to_compare = [xml_path_value]
            if 'extra_content_to_compare' in d.keys():
                content_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_content_to_compare'].items()]
            if 'extra_attrib_to_compare' in d.keys():
                attrib_to_compare += [XMLPathValue('PhysiCell_settings//' + k, v) for k, v in d['extra_attrib_to_compare'].items()]
            terminal_paths = []
            if 'terminal_paths' in d.keys():
                terminal_paths += d['terminal_paths']
            nonexistant_paths = []
            if 'nonexistant_paths' in d.keys():
                nonexistant_paths += d['nonexistant_paths']
            assert self.compare_xml(self.path_to_xml, "./config/__test__.xml", content_to_compare, attrib_to_compare, terminal_paths, nonexistant_paths, show_log=False)

    ########### perturb config tab in gui and check xml ############
    def create_config_tab_widget_dicts(self):
        self.create_config_tab_text_fields()
        self.create_config_tab_checkboxes()
        self.create_config_tab_radiobuttons()
        self.create_config_tab_comboboxes()

    def create_config_tab_text_fields(self):
        self.config_tab_text_fields = {'xmin': 'domain//x_min', 'xmax': 'domain//x_max', 'ymin': 'domain//y_min', 'ymax': 'domain//y_max', 'zmin': 'domain//z_min', 'zmax': 'domain//z_max', 'xdel': 'domain//dx', 'ydel': 'domain//dy', 'zdel': 'domain//dz',
                                'max_time': 'overall//max_time', 'diffusion_dt': 'overall//dt_diffusion', 'mechanics_dt': 'overall//dt_mechanics', 'phenotype_dt': 'overall//dt_phenotype',
                                'num_threads': 'parallel//omp_num_threads', 'folder': 'save//folder', 'full_interval': 'save//full_data//interval', 'svg_interval': 'save//SVG//interval',
                                'svg_substrate_min': 'save//SVG//plot_substrate//min_conc', 'svg_substrate_max': 'save//SVG//plot_substrate//max_conc', 'random_seed_integer': 'options//random_seed',
                                'csv_folder': 'initial_conditions//cell_positions//folder', 'csv_file': 'initial_conditions//cell_positions//filename'}
        self.config_tab_text_fields = {k: {'xml_path': v, 'new_value': "12321"} for k, v in self.config_tab_text_fields.items()}
        self.config_tab_text_fields["zmin"]["new_value"] = "0.12321"
        self.config_tab_text_fields["zmax"]["new_value"] = "0.12321"
        self.config_tab_text_fields["full_interval"]["extra_content_to_compare"] = {"save//SVG//interval": "12321"}
        self.config_tab_text_fields["svg_interval"]["extra_content_to_compare"] = {"save//full_data//interval": "12321"}
        if not self.xml_creator.config_tab.cells_csv.isChecked():
            pre_set_fn = lambda: self.xml_creator.config_tab.cells_csv.setChecked(True)
            post_set_fn = lambda: self.xml_creator.config_tab.cells_csv.setChecked(False)
            for field in ["csv_folder", "csv_file"]:
                self.config_tab_text_fields[field]["pre_set_fn"] = pre_set_fn
                self.config_tab_text_fields[field]["extra_attrib_to_compare"] = {"initial_conditions//cell_positions:enabled": "true"}
                self.config_tab_text_fields[field]["post_check_fn"] = post_set_fn
        if not self.xml_creator.config_tab.plot_substrate_svg.isChecked():
            new_extra_content_dict = {"save//SVG//plot_substrate//colormap": "YlOrRd", "save//SVG//plot_substrate//max_conc": None}
            if "extra_content_to_compare" not in self.config_tab_text_fields["svg_substrate_min"].keys():
                self.config_tab_text_fields["svg_substrate_min"]["extra_content_to_compare"] = new_extra_content_dict
            else:
                self.config_tab_text_fields["svg_substrate_min"]["extra_content_to_compare"].update(new_extra_content_dict)
            new_extra_content_dict = {"save//SVG//plot_substrate//colormap": "YlOrRd", "save//SVG//plot_substrate//min_conc": None}
            if "extra_content_to_compare" not in self.config_tab_text_fields["svg_substrate_max"].keys():
                self.config_tab_text_fields["svg_substrate_max"]["extra_content_to_compare"] = new_extra_content_dict
            else:
                self.config_tab_text_fields["svg_substrate_max"]["extra_content_to_compare"].update(new_extra_content_dict)
            limits_already_enabled = self.xml_creator.config_tab.plot_substrate_limits.isChecked()
            if limits_already_enabled:
                pre_set_fn = lambda: self.xml_creator.config_tab.plot_substrate_svg.setChecked(True)
                post_set_fn = lambda: self.xml_creator.config_tab.plot_substrate_svg.setChecked(False)
            else:
                pre_set_fn = lambda: (self.xml_creator.config_tab.plot_substrate_svg.setChecked(True), self.xml_creator.config_tab.plot_substrate_limits.setChecked(True))
                post_set_fn = lambda: (self.xml_creator.config_tab.plot_substrate_svg.setChecked(False), self.xml_creator.config_tab.plot_substrate_limits.setChecked(False))
            for field in ["svg_substrate_min", "svg_substrate_max"]:
                self.config_tab_text_fields[field]["pre_set_fn"] = pre_set_fn
                self.config_tab_text_fields[field]["extra_attrib_to_compare"] = {"save//SVG//plot_substrate:enabled": "true"}
                self.config_tab_text_fields[field]["post_check_fn"] = post_set_fn

    def create_config_tab_checkboxes(self):
        self.config_tab_checkboxes = {'save_full': 'save//full_data//enable', 'save_svg': 'save//SVG//enable', 
                             'plot_substrate_svg': 'save//SVG//plot_substrate:enabled', 
                             'plot_substrate_limits': 'save//SVG//plot_substrate:limits', 
                             'virtual_walls': 'options//virtual_wall_at_domain_edge', 'cells_csv': 'initial_conditions//cell_positions:enabled'}
        self.config_tab_checkboxes = {k: {'xml_path': v, 'new_value': lambda original_value: str(not original_value).lower()} for k, v in self.config_tab_checkboxes.items()}
    
    def create_config_tab_radiobuttons(self):
        self.config_tab_radiobuttons = {'random_seed_integer_button': 'options//random_seed', 'random_seed_random_button': 'options//random_seed'}
        self.config_tab_radiobuttons = {k: {'xml_path': v} for k, v in self.config_tab_radiobuttons.items()}
        self.config_tab_radiobuttons["random_seed_integer_button"]["new_value"] = lambda original_value: "system_clock" if original_value else "12321"
        self.config_tab_radiobuttons["random_seed_random_button"]["new_value"] = lambda original_value: "12321" if original_value else "system_clock"
        self.config_tab_radiobuttons["random_seed_integer_button"]["toggle_fn"] = self.toggle_random_seed_gp
        self.config_tab_radiobuttons["random_seed_random_button"]["toggle_fn"] = self.toggle_random_seed_gp

    def create_config_tab_comboboxes(self):
        self.config_tab_comboboxes = {'svg_substrate_colormap_dropdown': 'save//SVG//plot_substrate//colormap'}
        self.config_tab_comboboxes = {k: {'xml_path': v} for k, v in self.config_tab_comboboxes.items()}
        if not self.xml_creator.config_tab.plot_substrate_svg.isChecked():
            self.config_tab_comboboxes['svg_substrate_colormap_dropdown']['pre_loop_fn'] = lambda: self.xml_creator.config_tab.plot_substrate_svg.setChecked(True)
            self.config_tab_comboboxes['svg_substrate_colormap_dropdown']['extra_attrib_to_compare'] = {"save//SVG//plot_substrate:enabled": "true"}
            self.config_tab_comboboxes['svg_substrate_colormap_dropdown']['post_loop_fn'] = lambda: self.xml_creator.config_tab.plot_substrate_svg.setChecked(False)

    def toggle_random_seed_gp(self):
        current_id = self.xml_creator.config_tab.random_seed_gp.checkedId()
        self.xml_creator.config_tab.random_seed_gp.button((current_id+1) % 2).click()

    ########### perturb celldef tab in gui and check xml ############
    def create_celldef_tab_widget_dicts(self):
        # self.create_celldef_tab_text_fields()
        # self.create_celldef_tab_checkboxes()
        # self.create_celldef_tab_comboboxes()
        # self.create_celldef_tab_radiobuttons()
        self.create_celldef_tab_trees()
        self.create_celldef_tab_buttons()

    def create_celldef_tab_text_fields(self):
        raise NotImplementedError
    
    def create_celldef_tab_checkboxes(self):
        raise NotImplementedError
    
    def create_celldef_tab_comboboxes(self):
        raise NotImplementedError
    
    def create_celldef_tab_radiobuttons(self):
        raise NotImplementedError
    
    def create_celldef_tab_trees(self):
        self.celldef_tab_trees = {'tree': 'cell_definitions//cell_definition:name'}
        self.celldef_tab_trees = {k: {'xml_path': v} for k, v in self.celldef_tab_trees.items()}
    
    def create_celldef_tab_buttons(self):
        self.celldef_tab_buttons = {}
        new_button_new_value = "new_cell"
        self.celldef_tab_buttons['new_button'] = {}
        self.celldef_tab_buttons['new_button']['new_value'] = {} # no new_values to check. just stop recursing on getting to the new cell def
        self.celldef_tab_buttons['new_button']['post_push_fn'] = lambda: self.set_current_celldef_name(new_button_new_value)
        self.celldef_tab_buttons['new_button']['terminal_paths'] = [XMLPathValue(f'PhysiCell_settings//cell_definitions//cell_definition:name:{new_button_new_value}', None)]

        self.celldef_tab_buttons['delete_button'] = {'nonexistant_paths': [XMLPathValue(f'PhysiCell_settings//cell_definitions//cell_definition:name:{new_button_new_value}', None)]}

    def set_current_celldef_name(self, name):
        with suppress_output():
            self.xml_creator.celldef_tab.tree.currentItem().setText(0, name)

    def __del__(self):
        self.studio_app.quit()

if __name__ == "__main__":
    with suppress_output():
        studio_test = StudioTest(config_file="./config/PhysiCell_settings.xml")
    studio_test.check_import()
    # search through this file for all 
    studio_test.test_tab("config_tab", ["config_tab"])
    studio_test.test_tab("celldef_tab", ["celldef_tab"])
    