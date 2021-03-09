<template>
  <b-tab :id="scope+'-detector'" title="Detector">
    <div class="column">
      <fieldset>
        <legend>General</legend>
        <InputElement :textarea="data.name"
                      :id="scope+'-detector-name'"
                      title="Name" />
        <InputElement :time="data.modification_time"
                      :id="scope+'-detector-modification_time'"
                      title="Modified" />
        <InputElement :textarea="data.description"
                      :id="scope+'-detector-description'"
                      title="Description" />
        <InputElement :options="['TOF']"
                      :selected="data.detector_type"
                      :id="scope+'-detector-detector_type'"
                      title="Type" />
        <InputElement :input="data.timeres"
                      :id="scope+'-detector-timeres'"
                      title="Time res [ps]"
                      regex="double" />
        <div class="settings-row">
          <label class="settings-label">Virtual size:</label>
          <label>x:
            <input :id="scope+'-detector-0'" v-on:keypress="checkInput($event, 'double')"
                   class="settings-input" type="number" step="0.001" :value="data.virtual_size[0]">
          </label>
          <label>y:
          <input :id="scope+'-detector-1'" v-on:keypress="checkInput($event, 'double')"
                 class="settings-input" type="number" step="0.001" :value="data.virtual_size[1]">
          </label>
        </div>
      </fieldset>
    </div>
    <div class="column">
      <fieldset>
        <legend>Calibration settings</legend>
        <InputElement :input="data.tof_slope"
                      :id="scope+'-detector-tof_slope'"
                      title="ToF slope [s/channel]"
                      regex="scientific" />
        <InputElement :input="data.tof_offset"
                      :id="scope+'-detector-tof_offset'"
                      title="ToF offset [s]"
                      regex="scientific" />
        <InputElement :input="data.angle_slope"
                      :id="scope+'-detector-angle_slope'"
                      title="Angle slope [rad/channel]"
                      regex="double" />
        <InputElement :input="data.angle_offset"
                      :id="scope+'-detector-angle_offset'"
                      title="Angle offset [rad]"
                      regex="double" />
      </fieldset>
      <fieldset title="Efficiencies">
        <legend>Add efficiency files</legend>
      </fieldset>
    </div>
    <fieldset>
      <p>N.B Timing foils need to be carbon equivalent in numbers and have only one layer
        in order to calculate measurement energy spectra and depth profiles correctly.
      </p>
      <div class="detector-row">
        <span>Name of foil</span>
        <span>Distance from target</span>
        <span>Distance from previous foil</span>
        <span>Type</span>
        <span>Size</span>
        <span>Transmission</span>
      </div>
      <Foil v-for="foil in data.foils"
            :key="foil.name"
            :foil="foil"/>
    </fieldset>
  </b-tab>
</template>

<script>
  import InputElement from '../Elements/InputElement.vue';
  import Foil from '../Elements/Foil.vue';
  import Utils from '../../../js/utils';

  export default {
    name: 'Detector',
    components: {
      Foil,
      InputElement,
    },
    props: [
      'data',
      'scope',
    ],
    methods: {
      checkInput(keypress, regex) {
        Utils.checkInput(keypress, regex);
      },
    },
  };
</script>

<style>
  #Request-detector, #Measurement-detector {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    flex-wrap: wrap;
  }

  .detector-row {
    display: flex;
    justify-content: space-between;
    margin-bottom: 1em;
  }

  .detector-row input, span {
    width: 15%;
  }
</style>
