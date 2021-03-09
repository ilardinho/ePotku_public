<template>
  <b-tab :id="scope+'-mcsimu'" title="Simulation">
    <div class="column">
      <fieldset>
        <legend>General</legend>
        <InputElement :textarea="data.name"
                      :id="scope+'-mcsimu-name'"
                      title="Name" />
        <InputElement :time="data.modification_time"
                      :id="scope+'-mcsimu-modification_time'"
                      title="Modified" />
        <InputElement :textarea="data.description"
                      :id="scope+'-mcsimu-description'"
                      title="Description" />
      </fieldset>
      <fieldset>
        <legend>Physical parameters</legend>
        <InputElement :input="data.minimum_scattering_angle"
                      :id="scope+'-mcsimu-minimum_scattering_angle'"
                      title="Minimum scatter angle [°]"
                      regex="double" />
        <InputElement :input="data.minimum_main_scattering_angle"
                      :id="scope+'-mcsimu-minimum_main_scattering_angle'"
                      title="Minimum main scatter angle [°]"
                      regex="double" />
        <InputElement :input="data.minimum_energy"
                      :id="scope+'-mcsimu-minimum_energy'"
                      title="Minimum energy [MeV]"
                      regex="double" />
      </fieldset>
    </div>
    <div :id="scope+'-mcsimu-general-parameters'" class="column">
      <fieldset>
        <legend>General parameters</legend>
        <InputElement :options="['ERD', 'RBS']"
                      :selected="data.simulation_type"
                      :id="scope+'-mcsimu-simulation_type'"
                      title="Type of simulation"/>
        <!-- TODO: REC and SCT are saved as 'ERD' and 'RBS' in the settings files??? -->
        <InputElement :options="['Narrow', 'Wide']"
                      :selected="data.simulation_mode"
                      :id="scope+'-mcsimu-simulation_mode'"
                      title="Mode"/>
        <InputElement :input="data.number_of_ions"
                      :id="scope+'-mcsimu-number_of_ions'"
                      title="Number of ions"
                      regex="int" />
        <InputElement :input="data.number_of_preions"
                      :id="scope+'-mcsimu-number_of_preions'"
                      title="Number of pre-simulation ions"
                      regex="int" />
        <InputElement :input="data.seed_number"
                      :id="scope+'-mcsimu-seed_number'"
                      title="Seed for random number generator"
                      regex="int" />
        <InputElement :input="data.number_of_recoils"
                      :id="scope+'-mcsimu-number_of_recoils'"
                      title="Average number of recoils per primary ion"
                      regex="int" />
        <InputElement :input="data.number_of_scaling_ions"
                      :id="scope+'-mcsimu-number_of_scaling_ions'"
                      title="Number of real ions per each scaling ion"
                      regex="int" />
      </fieldset>
    </div>
  </b-tab>
</template>

<script>
  import InputElement from '../Elements/InputElement.vue';

  export default {
    name: 'Simulation',
    components: {
      InputElement,
    },
    props: [
      'data',
      'scope',
    ],
  };
</script>

<style>
  #Request-mcsimu, #Measurement-mcsimu {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
  }
</style>
