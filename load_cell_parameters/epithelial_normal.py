import xml.etree.ElementTree as ET

def get_epithelial_normal_template():
    epithelial_normal_template = '''
    <cell_definitions>
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
                    <cell_adhesion_affinities>
                        <cell_adhesion_affinity name="fibroblast">1.0</cell_adhesion_affinity>
                        <cell_adhesion_affinity name="epithelial_normal">1</cell_adhesion_affinity>
                        <cell_adhesion_affinity name="epithelial_tumor">1.0</cell_adhesion_affinity>
                        <cell_adhesion_affinity name="mesenchymal_normal">1.0</cell_adhesion_affinity>
                        <cell_adhesion_affinity name="mesenchymal_tumor">1.0</cell_adhesion_affinity>
                    </cell_adhesion_affinities>
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
                <secretion>
                    <substrate name="inflammatory_signal">
                        <secretion_rate units="1/min">0</secretion_rate>
                        <secretion_target units="substrate density">1</secretion_target>
                        <uptake_rate units="1/min">0</uptake_rate>
                        <net_export_rate units="total substrate/min">0</net_export_rate>
                    </substrate>
                    <substrate name="ecm">
                        <secretion_rate units="1/min">0.0</secretion_rate>
                        <secretion_target units="substrate density">0.0</secretion_target>
                        <uptake_rate units="1/min">0.0</uptake_rate>
                        <net_export_rate units="total substrate/min">0.0</net_export_rate>
                    </substrate>
                </secretion>
                <cell_interactions>
                    <dead_phagocytosis_rate units="1/min">0</dead_phagocytosis_rate>
                    <live_phagocytosis_rates>
                        <phagocytosis_rate name="fibroblast" units="1/min">0.0</phagocytosis_rate>
                        <phagocytosis_rate name="epithelial_normal" units="1/min">0</phagocytosis_rate>
                        <phagocytosis_rate name="epithelial_tumor" units="1/min">0.0</phagocytosis_rate>
                        <phagocytosis_rate name="mesenchymal_normal" units="1/min">0.0</phagocytosis_rate>
                        <phagocytosis_rate name="mesenchymal_tumor" units="1/min">0.0</phagocytosis_rate>
                    </live_phagocytosis_rates>
                    <attack_rates>
                        <attack_rate name="fibroblast" units="1/min">0.0</attack_rate>
                        <attack_rate name="epithelial_normal" units="1/min">0</attack_rate>
                        <attack_rate name="epithelial_tumor" units="1/min">0.0</attack_rate>
                        <attack_rate name="mesenchymal_normal" units="1/min">0.0</attack_rate>
                        <attack_rate name="mesenchymal_tumor" units="1/min">0.0</attack_rate>
                    </attack_rates>
                    <damage_rate units="1/min">1</damage_rate>
                    <fusion_rates>
                        <fusion_rate name="fibroblast" units="1/min">0.0</fusion_rate>
                        <fusion_rate name="epithelial_normal" units="1/min">0</fusion_rate>
                        <fusion_rate name="epithelial_tumor" units="1/min">0.0</fusion_rate>
                        <fusion_rate name="mesenchymal_normal" units="1/min">0.0</fusion_rate>
                        <fusion_rate name="mesenchymal_tumor" units="1/min">0.0</fusion_rate>
                    </fusion_rates>
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates>
                        <transformation_rate name="epithelial_normal" units="1/min">0</transformation_rate>
                        <transformation_rate name="fibroblast" units="1/min">0.0</transformation_rate>
                        <transformation_rate name="epithelial_tumor" units="1/min">0.0</transformation_rate>
                        <transformation_rate name="mesenchymal_normal" units="1/min">0.0</transformation_rate>
                        <transformation_rate name="mesenchymal_tumor" units="1/min">0.0</transformation_rate>
                    </transformation_rates>
                </cell_transformations>
            </phenotype>
            <custom_data>
                <sample conserved="false" units="dimensionless" description="">1.0</sample>
            </custom_data>
    '''

    return epithelial_normal_template


    # #takes string containing XML data and converts into an Element object- returns the root element of the XML tree s.t. access child elements
    # root = ET.fromstring(epithelial_normal_template)
    # return root
