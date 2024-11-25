
import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html
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
                        <chemotaxis>
                            <enabled>false</enabled>
                            <substrate>inflammatory_signal</substrate>
                            <direction>1</direction>
                        </chemotaxis>
                        <advanced_chemotaxis>
                            <enabled>false</enabled>
                            <normalize_each_gradient>false</normalize_each_gradient>
                            <chemotactic_sensitivities>
                                <chemotactic_sensitivity substrate="inflammatory_signal">0.0</chemotactic_sensitivity>
                                <chemotactic_sensitivity substrate="ecm">0.0</chemotactic_sensitivity>
                            </chemotactic_sensitivities>
                        </advanced_chemotaxis>
                    </options>
                </motility>
                <secretion />
                <cell_interactions>
                    <apoptotic_phagocytosis_rate units="1/min">0</apoptotic_phagocytosis_rate>
                    <necrotic_phagocytosis_rate units="1/min">0</necrotic_phagocytosis_rate>
                    <other_dead_phagocytosis_rate units="1/min">0</other_dead_phagocytosis_rate>
                    <live_phagocytosis_rates />
                    <attack_rates />
                    <damage_rate units="1/min">1</damage_rate>
                    <fusion_rates />
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates />
                </cell_transformations>
            </phenotype>
    '''
x = ET.fromstring(epithelial_normal_template)

