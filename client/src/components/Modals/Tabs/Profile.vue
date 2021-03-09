<template>
  <b-tab :id="scope+'-profile'" title="Profile">
    <div class="column">
      <fieldset>
        <legend>General</legend>
        <InputElement :textarea="data.general.name"
                      :id="scope+'-profile-name'"
                      title="Name" />
        <InputElement :time="data.general.modification_time"
                      :id="scope+'-profile-modification_time'"
                      title="Modified" />
        <InputElement :textarea="data.general.description"
                      :id="scope+'-profile-description'"
                      title="Description" />
      </fieldset>
      <fieldset>
        <legend>Energy spectra</legend>
        <InputElement :input="data.energy_spectra.channel_width"
                      :id="scope+'-profile-channel_width'"
                      title="Channel width [MeV]"
                      regex="double"/>
      </fieldset>
      <fieldset>
        <legend>Composition changes</legend>
        <InputElement :input="data.composition_changes.number_of_splits"
                      :id="scope+'-profile-number_of_splits'"
                      title="Number of splits"
                      regex="int"/>
        <InputElement :options="['First']"
                      :selected="data.composition_changes.normalization"
                      :id="scope+'-profile-normalization'"
                      title="Normalization" />
      </fieldset>
    </div>
    <div class="column">
      <fieldset>
        <legend>Depth profile</legend>
        <InputElement :input="data.depth_profiles.reference_density"
                      :id="scope+'-profile-reference_density'"
                      title="Reference density [g/cm3]"
                      regex="double"/>
        <InputElement :input="data.depth_profiles.number_of_depth_steps"
                      :id="scope+'-profile-number_of_depth_steps'"
                      title="Number of depth steps"
                      regex="int" />
        <InputElement :input="data.depth_profiles.depth_step_for_stopping"
                      :id="scope+'-profile-depth_step_for_stopping'"
                      title="Step for stopping [1e15 at./cm2]"
                      regex="int" />
        <InputElement :input="data.depth_profiles.depth_step_for_output"
                      :id="scope+'-profile-depth_step_for_output'"
                      title="Step for stopping [1e15 at./cm2]"
                      regex="int" />
        <div class="settings-row">
          <label class="settings-label">Concentration scaling 100% [1e15 at./cm2]:</label>
          <label class="settings-label">From:
            <input :id="scope+'-profile-depth_for_concentration_from'" class="settings-input"
                   :value="data.depth_profiles.depth_for_concentration_from"
                   v-on:keypress="checkInput($event, 'int')"
                   type="number">
          </label>
          <label class="settings-label">To:
            <input :id="scope+'-profile-depth_for_concentration_to'" class="settings-input"
                   :value="data.depth_profiles.depth_for_concentration_to"
                   v-on:keypress="checkInput($event, 'int')"
                   type="number">
          </label>
        </div>
      </fieldset>
    </div>
  </b-tab>
</template>

<script>
  import InputElement from '../Elements/InputElement.vue';
  import Utils from '../../../js/utils';

  export default {
    name: 'Profile',
    props: [
      'data',
      'scope',
    ],
    components: {
      InputElement,
    },
    methods: {
      checkInput(keypress, regex) {
        Utils.checkInput(keypress, regex);
      },
    },
  };
</script>

<style>
  #Request-profile, #Measurement-profile{
    display: flex;
    flex-direction: row;
    justify-content: space-between;
  }
</style>
