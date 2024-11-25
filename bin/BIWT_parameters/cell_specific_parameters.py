def get_default():
    default_template = '''
        <phenotype>
                    <cycle code="6" name="Flow cytometry model (separated)">
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="false">300.</duration>
                            <duration index="1" fixed_duration="true">480.</duration>
                            <duration index="2" fixed_duration="true">240.</duration>
                            <duration index="3" fixed_duration="true">60.</duration>
                        </phase_durations>
                    </cycle>
                    <death>
                        <model code="100" name="apoptosis">
                            <death_rate units="1/min">5.31667e-05</death_rate>
                            <phase_durations units="min">
                                <duration index="0" fixed_duration="true">516</duration>
                            </phase_durations>
                            <parameters>
                                <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                                <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                                <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                                <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                                <calcification_rate units="1/min">0</calcification_rate>
                                <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                            </parameters>
                        </model>
                        <model code="101" name="necrosis">
                            <death_rate units="1/min">0.0</death_rate>
                            <phase_durations units="min">
                                <duration index="0" fixed_duration="true">0</duration>
                                <duration index="1" fixed_duration="true">86400</duration>
                            </phase_durations>
                            <parameters>
                                <unlysed_fluid_change_rate units="1/min">1.11667e-2</unlysed_fluid_change_rate>
                                <lysed_fluid_change_rate units="1/min">8.33333e-4</lysed_fluid_change_rate>
                                <cytoplasmic_biomass_change_rate units="1/min">5.33333e-5</cytoplasmic_biomass_change_rate>
                                <nuclear_biomass_change_rate units="1/min">2.16667e-3</nuclear_biomass_change_rate>
                                <calcification_rate units="1/min">0</calcification_rate>
                                <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                            </parameters>
                        </model>
                    </death>
                    <volume>
                        <total units="micron^3">2494</total>
                        <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                        <nuclear units="micron^3">540</nuclear>
                        <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                        <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                        <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                        <calcified_fraction units="dimensionless">0</calcified_fraction>
                        <calcification_rate units="1/min">0</calcification_rate>
                        <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                    </volume>
                    <mechanics>
                        <cell_cell_adhesion_strength units="micron/min">0.4</cell_cell_adhesion_strength>
                        <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                        <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
                        <cell_adhesion_affinities>
                            <cell_adhesion_affinity name="default">1</cell_adhesion_affinity>
                        </cell_adhesion_affinities>
                        <options>
                            <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                            <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                        </options>
                        <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                        <attachment_rate units="1/min">0.0</attachment_rate>
                        <detachment_rate units="1/min">0.0</detachment_rate>
                        <maximum_number_of_attachments>12</maximum_number_of_attachments>
                    </mechanics>
                    <motility>
                        <speed units="micron/min">1</speed>
                        <persistence_time units="min">1</persistence_time>
                        <migration_bias units="dimensionless">.5</migration_bias>
                        <options>
                            <enabled>false</enabled>
                            <use_2D>true</use_2D>                    
                        </options>
                    </motility>
                    <secretion>
                        <substrate name="substrate">
                            <secretion_rate units="1/min">0</secretion_rate>
                            <secretion_target units="substrate density">1</secretion_target>
                            <uptake_rate units="1/min">0</uptake_rate>
                            <net_export_rate units="total substrate/min">0</net_export_rate>
                        </substrate>
                    </secretion>
                    <cell_interactions>
                        <apoptotic_phagocytosis_rate units="1/min">0.0</apoptotic_phagocytosis_rate>
                        <necrotic_phagocytosis_rate units="1/min">0.0</necrotic_phagocytosis_rate>
                        <other_dead_phagocytosis_rate units="1/min">0.0</other_dead_phagocytosis_rate>
                        <live_phagocytosis_rates>
                            <phagocytosis_rate name="default" units="1/min">0</phagocytosis_rate>
                        </live_phagocytosis_rates>
                        <attack_rates>
                            <attack_rate name="default" units="1/min">0</attack_rate>
                        </attack_rates>
                        <attack_damage_rate units="1/min">1</attack_damage_rate>
                        <attack_duration units="min">0.1</attack_duration>
                        <fusion_rates>
                            <fusion_rate name="default" units="1/min">0</fusion_rate>
                        </fusion_rates>
                    </cell_interactions>
                    <cell_transformations>
                        <transformation_rates>
                            <transformation_rate name="default" units="1/min">0</transformation_rate>
                        </transformation_rates>
                    </cell_transformations>
                    <cell_integrity>
                        <damage_rate units="1/min">0.0</damage_rate>
                        <damage_repair_rate units="1/min">0.0</damage_repair_rate>
                    </cell_integrity>
                </phenotype>
        '''
    return default_template

# start of epi-caf-hypothesis1
def get_epithelial_normal():
    epithelial_normal_template = '''
            <phenotype>
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">0</rate>
                    </phase_transition_rates>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">5.31667e-05</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">516</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0.0</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">0</duration>
                            <duration index="1" fixed_duration="true">86400</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">1.11667e-2</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">8.33333e-4</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">5.33333e-5</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">2.16667e-3</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <total units="micron^3">2494</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0</calcified_fraction>
                    <calcification_rate units="1/min">0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">0.4</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
                    <cell_adhesion_affinities />
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0.0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                </mechanics>
                <motility>
                    <speed units="micron/min">0</speed>
                    <persistence_time units="min">1</persistence_time>
                    <migration_bias units="dimensionless">.5</migration_bias>
                    <options>
                        <enabled>false</enabled>
                        <use_2D>true</use_2D>
                    </options>
                </motility>
                <secretion />
                <cell_interactions>
                    <apoptotic_phagocytosis_rate units="1/min">0</apoptotic_phagocytosis_rate>
                    <necrotic_phagocytosis_rate units="1/min">0</necrotic_phagocytosis_rate>
                    <other_dead_phagocytosis_rate units="1/min">0</other_dead_phagocytosis_rate>
                    <live_phagocytosis_rates />
                    <attack_rates />
                    <attack_damage_rate units="1/min">1</attack_damage_rate>
                    <fusion_rates />
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />
                </cell_transformations>
            </phenotype>
        '''

    return epithelial_normal_template

def get_mesenchymal_normal():
    mesenchymal_normal_template = '''
        <phenotype>
                    <cycle code="5" name="live">
                        <phase_transition_rates units="1/min">
                            <rate start_index="0" end_index="0" fixed_duration="false">0</rate>
                        </phase_transition_rates>
                    </cycle>
                    <death>
                        <model code="100" name="apoptosis">
                            <death_rate units="1/min">5.31667e-06</death_rate>
                            <phase_durations units="min">
                                <duration index="0" fixed_duration="true">516</duration>
                            </phase_durations>
                            <parameters>
                                <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                                <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                                <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                                <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                                <calcification_rate units="1/min">0</calcification_rate>
                                <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                            </parameters>
                        </model>
                        <model code="101" name="necrosis">
                            <death_rate units="1/min">0.0</death_rate>
                            <phase_durations units="min">
                                <duration index="0" fixed_duration="true">0</duration>
                                <duration index="1" fixed_duration="true">86400</duration>
                            </phase_durations>
                            <parameters>
                                <unlysed_fluid_change_rate units="1/min">1.11667e-2</unlysed_fluid_change_rate>
                                <lysed_fluid_change_rate units="1/min">8.33333e-4</lysed_fluid_change_rate>
                                <cytoplasmic_biomass_change_rate units="1/min">5.33333e-5</cytoplasmic_biomass_change_rate>
                                <nuclear_biomass_change_rate units="1/min">2.16667e-3</nuclear_biomass_change_rate>
                                <calcification_rate units="1/min">0</calcification_rate>
                                <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                            </parameters>
                        </model>
                    </death>
                    <volume>
                        <total units="micron^3">2494</total>
                        <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                        <nuclear units="micron^3">540</nuclear>
                        <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                        <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                        <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                        <calcified_fraction units="dimensionless">0</calcified_fraction>
                        <calcification_rate units="1/min">0</calcification_rate>
                        <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                    </volume>
                    <mechanics>
                        <cell_cell_adhesion_strength units="micron/min">0.2</cell_cell_adhesion_strength>
                        <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                        <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
                        <cell_adhesion_affinities />                           
                        <options>
                            <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                            <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                        </options>
                        <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                        <attachment_rate units="1/min">0.0</attachment_rate>
                        <detachment_rate units="1/min">0.0</detachment_rate>
                    </mechanics>
                    <motility>
                        <speed units="micron/min">0.249908679</speed>
                        <persistence_time units="min">10</persistence_time>
                        <migration_bias units="dimensionless">.5</migration_bias>
                        <options>
                            <enabled>true</enabled>
                            <use_2D>true</use_2D>
                        </options>
                    </motility>
                    <secretion/>
                    <cell_interactions>
                        <dead_phagocytosis_rate units="1/min">0</dead_phagocytosis_rate>
                        <live_phagocytosis_rates />
                        <attack_rates />            
                        <attack_damage_rate units="1/min">1</attack_damage_rate>
                        <fusion_rates />                           
                    </cell_interactions>
                    <cell_transformations>
                        <transformation_rates />                            
                    </cell_transformations>
                </phenotype>
     '''
    return mesenchymal_normal_template

def get_fibroblast():
    fibroblast_template = '''
        <phenotype>
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">0</rate>
                    </phase_transition_rates>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">0</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">516</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0.0</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">0</duration>
                            <duration index="1" fixed_duration="true">86400</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">1.11667e-2</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">8.33333e-4</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">5.33333e-5</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">2.16667e-3</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <total units="micron^3">2494</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0</calcified_fraction>
                    <calcification_rate units="1/min">0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">0.4</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
                    <cell_adhesion_affinities/>
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0.0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                </mechanics>
                <motility>
                    <speed units="micron/min">8.03971290256818e-06</speed>
                    <persistence_time units="min">10</persistence_time>
                    <migration_bias units="dimensionless">.5</migration_bias>
                    <options>
                        <enabled>true</enabled>
                        <use_2D>true</use_2D>
                    </options>
                </motility>
                <secretion />                           
                <cell_interactions>
                    <dead_phagocytosis_rate units="1/min">0</dead_phagocytosis_rate>
                    <live_phagocytosis_rates />
                    <attack_rates />
                    <attack_damage_rate units="1/min">1</attack_damage_rate>
                    <fusion_rates />
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />
                </cell_transformations>
            </phenotype>
    '''
    return fibroblast_template

def get_tumor_epithlial():
    epithilial_tumor_template = '''
        <phenotype>
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">1e-3</rate>
                    </phase_transition_rates>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">5.31667e-05</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">516</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0.0</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">0</duration>
                            <duration index="1" fixed_duration="true">86400</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">1.11667e-2</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">8.33333e-4</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">5.33333e-5</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">2.16667e-3</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <total units="micron^3">2494</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0</calcified_fraction>
                    <calcification_rate units="1/min">0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">0.4</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
                    <cell_adhesion_affinities />                    
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0.0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                </mechanics>
                <motility>
                    <speed units="micron/min">0</speed>
                    <persistence_time units="min">1</persistence_time>
                    <migration_bias units="dimensionless">.5</migration_bias>
                    <options>
                        <enabled>false</enabled>
                        <use_2D>true</use_2D>
                    </options>
                </motility>
                <secretion />
                <cell_interactions>
                    <dead_phagocytosis_rate units="1/min">0</dead_phagocytosis_rate>
                    <live_phagocytosis_rates />
                    <attack_rates />                        
                    <attack_damage_rate units="1/min">1</attack_damage_rate>
                    <fusion_rates />                
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />
                </cell_transformations>
            </phenotype>
    '''
    return epithilial_tumor_template

def get_mesenchymal_tumor():
    mesenchymal_tumor_template = '''
        <phenotype>
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">0</rate>
                    </phase_transition_rates>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">5.31667e-06</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">516</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0.0</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">0</duration>
                            <duration index="1" fixed_duration="true">86400</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">1.11667e-2</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">8.33333e-4</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">5.33333e-5</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">2.16667e-3</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <total units="micron^3">2494</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0</calcified_fraction>
                    <calcification_rate units="1/min">0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">0.2</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
                    <cell_adhesion_affinities />
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0.0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                </mechanics>
                <motility>
                    <speed units="micron/min">0.249908679</speed>
                    <persistence_time units="min">10</persistence_time>
                    <migration_bias units="dimensionless">.5</migration_bias>
                    <options>
                        <enabled>true</enabled>
                        <use_2D>true</use_2D>
                    </options>
                </motility>
                <secretion />
                <cell_interactions>
                    <dead_phagocytosis_rate units="1/min">0</dead_phagocytosis_rate>
                    <live_phagocytosis_rates />
                    <attack_rates />                        
                    <attack_damage_rate units="1/min">1</attack_damage_rate>
                    <fusion_rates />                
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />
                </cell_transformations>
            </phenotype>
    '''
    return mesenchymal_tumor_template
#end of epi-caf-hypothesis1

#start of example3_immune, skipping tumor
def get_machrophage():
    machrophage_template = '''
     <phenotype>
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">0</rate>
                    </phase_transition_rates>
                    <standard_asymmetric_division enabled="False" />
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">0</death_rate>
                        <phase_transition_rates units="1/min">
                            <rate start_index="0" end_index="1" fixed_duration="false">0.001938</rate>
                        </phase_transition_rates>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0.0</death_rate>
                        <phase_transition_rates units="1/min">
                            <rate start_index="0" end_index="1" fixed_duration="false">9000000000.0</rate>
                            <rate start_index="1" end_index="2" fixed_duration="true">1.15741e-05</rate>
                        </phase_transition_rates>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">1.11667e-02</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">8.33333e-4</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">5.33333e-05</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">2.16667e-4</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">7e-05</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <total units="micron^3">2494</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0.0</calcified_fraction>
                    <calcification_rate units="1/min">0.0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">0</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
                    <cell_adhesion_affinities />
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0.0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                    <maximum_number_of_attachments>12</maximum_number_of_attachments>
                </mechanics>
                <motility>
                    <speed units="micron/min">1.0</speed>
                    <persistence_time units="min">5</persistence_time>
                    <migration_bias units="dimensionless">0.25</migration_bias>
                    <options>
                        <enabled>true</enabled>
                        <use_2D>true</use_2D>
                    </options>
                </motility>
                <secretion />
                <cell_interactions>
                    <apoptotic_phagocytosis_rate units="1/min">0.1</apoptotic_phagocytosis_rate>
                    <necrotic_phagocytosis_rate units="1/min">0.017</necrotic_phagocytosis_rate>
                    <other_dead_phagocytosis_rate units="1/min">0.0</other_dead_phagocytosis_rate>
                    <live_phagocytosis_rates />
                    <attack_rates />                        
                    <attack_damage_rate units="1/min">1</attack_damage_rate>
                    <fusion_rates />                
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />
                </cell_transformations>
            </phenotype> 
    '''
    return  machrophage_template

def get_CD8_T_cell():
    CD8_T_Cell_Template = '''
        <phenotype>
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">0</rate>
                    </phase_transition_rates>
                    <standard_asymmetric_division enabled="False">
                        <asymmetric_division_probability name="tumor" units="dimensionless">0</asymmetric_division_probability>
                        <asymmetric_division_probability name="macrophage" units="dimensionless">0</asymmetric_division_probability>
                        <asymmetric_division_probability name="CD8 T cell" units="dimensionless">1.0</asymmetric_division_probability>
                    </standard_asymmetric_division>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">0</death_rate>
                        <phase_transition_rates units="1/min">
                            <rate start_index="0" end_index="1" fixed_duration="false">0.001938</rate>
                        </phase_transition_rates>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0.0</death_rate>
                        <phase_transition_rates units="1/min">
                            <rate start_index="0" end_index="1" fixed_duration="false">9000000000.0</rate>
                            <rate start_index="1" end_index="2" fixed_duration="true">1.15741e-05</rate>
                        </phase_transition_rates>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">1.11667e-02</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">8.33333e-4</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">5.33333e-05</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">2.16667e-4</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">7e-05</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <total units="micron^3">2494</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0.0</calcified_fraction>
                    <calcification_rate units="1/min">0.0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">0</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
                    <cell_adhesion_affinities />
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0.0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                    <maximum_number_of_attachments>12</maximum_number_of_attachments>
                </mechanics>
                <motility>
                    <speed units="micron/min">1.0</speed>
                    <persistence_time units="min">5</persistence_time>
                    <migration_bias units="dimensionless">0.5</migration_bias>
                    <options>
                        <enabled>true</enabled>
                        <use_2D>true</use_2D>
                    </options>
                </motility>
                <secretion />                
                <cell_interactions>
                    <apoptotic_phagocytosis_rate units="1/min">0.0</apoptotic_phagocytosis_rate>
                    <necrotic_phagocytosis_rate units="1/min">0.0</necrotic_phagocytosis_rate>
                    <other_dead_phagocytosis_rate units="1/min">0.0</other_dead_phagocytosis_rate>
                    <live_phagocytosis_rates />                       
                    <attack_rates />                      
                    <attack_damage_rate units="1/min">1.0</attack_damage_rate>
                    <attack_duration units="min">15</attack_duration>
                    <fusion_rates />                       
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />    
                 </cell_transformations>                   
                <cell_integrity>
                    <damage_rate units="1/min">0</damage_rate>
                    <damage_repair_rate units="1/min">0.0</damage_repair_rate>
                </cell_integrity>
            </phenotype>
    '''
    return CD8_T_Cell_Template
#end of example3_immune, skipping tumor

#start of example3_tam_egf_go
def get_M0_macrophage():
    M0_macrophage_template = '''
        <phenotype>
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">7.2e-5</rate>
                    </phase_transition_rates>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">7.2e-5</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">516</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">0</duration>
                            <duration index="1" fixed_duration="true">86400</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <total units="micron^3">2494</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0</calcified_fraction>
                    <calcification_rate units="1/min">0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">0</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.5</relative_maximum_adhesion_distance>
                    <cell_adhesion_affinities>
                        <cell_adhesion_affinity name="tumor">1.0</cell_adhesion_affinity>
                        <cell_adhesion_affinity name="M0 macrophage">1.0</cell_adhesion_affinity>
                        <cell_adhesion_affinity name="M1 macrophage">1.0</cell_adhesion_affinity>
                        <cell_adhesion_affinity name="M2 macrophage">1.0</cell_adhesion_affinity>
                        <cell_adhesion_affinity name="Th2 CD4 T cell">1.0</cell_adhesion_affinity>
                    </cell_adhesion_affinities>
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                    <maximum_number_of_attachments>12</maximum_number_of_attachments>
                </mechanics>
                <motility>
                    <speed units="micron/min">0.5</speed>
                    <persistence_time units="min">5</persistence_time>
                    <migration_bias units="dimensionless">0.25</migration_bias>
                    <options>
                        <enabled>true</enabled>
                        <use_2D>true</use_2D>
                    </options>
                </motility>
                <secretion />
                <cell_interactions>
                    <apoptotic_phagocytosis_rate units="1/min">0.001</apoptotic_phagocytosis_rate>
                    <necrotic_phagocytosis_rate units="1/min">0.001</necrotic_phagocytosis_rate>
                    <other_dead_phagocytosis_rate units="1/min">0.001</other_dead_phagocytosis_rate>
                    <live_phagocytosis_rates />                      
                    <attack_rates />                       
                    <attack_damage_rate units="1/min">1</attack_damage_rate>
                    <attack_duration units="min">0.1</attack_duration>
                    <fusion_rates />                    
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />                      
                </cell_transformations>
                <cell_integrity>
                    <damage_rate units="1/min">0.0</damage_rate>
                    <damage_repair_rate units="1/min">0.0</damage_repair_rate>
                </cell_integrity>
            </phenotype>
    '''
    return  M0_macrophage_template

def get_M1_macrophage():
    M1_macrophage_template = '''
        <phenotype>
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">7.2e-5</rate>
                    </phase_transition_rates>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">7.2e-5</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">516</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">0</duration>
                            <duration index="1" fixed_duration="true">86400</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <total units="micron^3">2494</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0</calcified_fraction>
                    <calcification_rate units="1/min">0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">0</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.5</relative_maximum_adhesion_distance>
                    <cell_adhesion_affinities />                      
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                    <maximum_number_of_attachments>12</maximum_number_of_attachments>
                </mechanics>
                <motility>
                    <speed units="micron/min">0.5</speed>
                    <persistence_time units="min">5</persistence_time>
                    <migration_bias units="dimensionless">0.25</migration_bias>
                    <options>
                        <enabled>true</enabled>
                        <use_2D>true</use_2D>
                    </options>
                </motility>
                <secretion />
                <cell_interactions>
                    <apoptotic_phagocytosis_rate units="1/min">0.001</apoptotic_phagocytosis_rate>
                    <necrotic_phagocytosis_rate units="1/min">0.001</necrotic_phagocytosis_rate>
                    <other_dead_phagocytosis_rate units="1/min">0.001</other_dead_phagocytosis_rate>
                    <live_phagocytosis_rates />                      
                    <attack_rates />                      
                    <attack_damage_rate units="1/min">1</attack_damage_rate>
                    <attack_duration units="min">0.1</attack_duration>
                    <fusion_rates />
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />
                </cell_transformations>
                <cell_integrity>
                    <damage_rate units="1/min">0.0</damage_rate>
                    <damage_repair_rate units="1/min">0.0</damage_repair_rate>
                </cell_integrity>
            </phenotype>
    '''
    return M1_macrophage_template

def get_M2_macrophage():
    M2_macrophage_template = '''
        <phenotype>
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">7.2e-5</rate>
                    </phase_transition_rates>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">7.2e-5</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">516</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">0</duration>
                            <duration index="1" fixed_duration="true">86400</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <total units="micron^3">2494</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0</calcified_fraction>
                    <calcification_rate units="1/min">0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">0</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.5</relative_maximum_adhesion_distance>
                    <cell_adhesion_affinities />       
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                    <maximum_number_of_attachments>12</maximum_number_of_attachments>
                </mechanics>
                <motility>
                    <speed units="micron/min">0.5</speed>
                    <persistence_time units="min">5</persistence_time>
                    <migration_bias units="dimensionless">0.25</migration_bias>
                    <options>
                        <enabled>true</enabled>
                        <use_2D>true</use_2D>
                    </options>
                </motility>
                <secretion />
                </secretion>
                <cell_interactions>
                    <apoptotic_phagocytosis_rate units="1/min">0.001</apoptotic_phagocytosis_rate>
                    <necrotic_phagocytosis_rate units="1/min">0.001</necrotic_phagocytosis_rate>
                    <other_dead_phagocytosis_rate units="1/min">0.001</other_dead_phagocytosis_rate>
                    <live_phagocytosis_rates />                       
                    <attack_rates />                        
                    <attack_damage_rate units="1/min">1</attack_damage_rate>
                    <attack_duration units="min">0.1</attack_duration>
                    <fusion_rates />                       
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />                      
                </cell_transformations>
                <cell_integrity>
                    <damage_rate units="1/min">0.0</damage_rate>
                    <damage_repair_rate units="1/min">0.0</damage_repair_rate>
                </cell_integrity>
            </phenotype>
    '''
    return M2_macrophage_template

def get_Th2_CD4_T_cell():
    Th2_CD4_T_cell_tempate = '''
        <phenotype>
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">0</rate>
                    </phase_transition_rates>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">7.2e-5</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">516</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0</death_rate>
                        <phase_durations units="min">
                            <duration index="0" fixed_duration="true">0</duration>
                            <duration index="1" fixed_duration="true">86400</duration>
                        </phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <total units="micron^3">2494</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0</calcified_fraction>
                    <calcification_rate units="1/min">0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">0</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.5</relative_maximum_adhesion_distance>
                    <cell_adhesion_affinities />                      
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                    <maximum_number_of_attachments>12</maximum_number_of_attachments>
                </mechanics>
                <motility>
                    <speed units="micron/min">0.5</speed>
                    <persistence_time units="min">5</persistence_time>
                    <migration_bias units="dimensionless">0.5</migration_bias>
                    <options>
                        <enabled>true</enabled>
                        <use_2D>true</use_2D>
                    </options>
                </motility>
                <secretion />                  
                <cell_interactions>
                    <apoptotic_phagocytosis_rate units="1/min">0.0</apoptotic_phagocytosis_rate>
                    <necrotic_phagocytosis_rate units="1/min">0.0</necrotic_phagocytosis_rate>
                    <other_dead_phagocytosis_rate units="1/min">0.0</other_dead_phagocytosis_rate>
                    <live_phagocytosis_rates />                     
                    <attack_rates />                      
                    <attack_damage_rate units="1/min">1</attack_damage_rate>
                    <attack_duration units="min">0.1</attack_duration>
                    <fusion_rates />                      
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />                      
                </cell_transformations>
                <cell_integrity>
                    <damage_rate units="1/min">0.0</damage_rate>
                    <damage_repair_rate units="1/min">0.0</damage_repair_rate>
                </cell_integrity>
            </phenotype>

    '''
    return Th2_CD4_T_cell_tempate
# end of example3_tam_egf_go


