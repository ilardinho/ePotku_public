<template>
  <b-tab :id="scope+'-measurement'" title="Measurement" active>
    <div class="column">
      <fieldset>
        <legend>General</legend>
        <InputElement :textarea="data.general.name"
                      :id="scope+'-measurement-name'"
                      title="Name" />
        <InputElement :time="data.general.modification_time"
                      :id="scope+'-measurement-modification_time'"
                      title="Modified" />
        <InputElement :textarea="data.general.description"
                      :id="scope+'-measurement-description'"
                      title="Description" />
      </fieldset>
      <fieldset>
        <legend>Measurement setup</legend>
        <img title="test" id="setup-angle" src="../../../assets/measurement_setup_angles.png">
      </fieldset>
      <fieldset>
        <legend>Geometry</legend>
        <InputElement :input="data.geometry.detector_theta"
                      :id="scope+'-measurement-detector_theta'"
                      title="Detector angle (θ) [°]"
                      regex="double" />
        <InputElement :input="0.0"
                      title="Detector angle (Φ) [°]"
                      :disabled="true" />
        <InputElement :input="data.geometry.target_theta"
                      :id="scope+'-measurement-target_theta'"
                      title="Target angle (β) [°]"
                      regex="double" />
        <InputElement :input="0.0"
                      title="Detector angle (Φ) [°]"
                      :disabled="true" />
      </fieldset>
    </div>
    <div class="column">
      <fieldset>
        <legend>Beam settings</legend>
        <InputElement :button="data.beam.element"
                      :id="scope+'-measurement-element'"
                      title="Beam ion"
                      regex="int" />
        <!-- TODO: dropdown list for valid isotopes -->
        <InputElement :input="data.beam.isotope"
                      :id="scope+'-measurement-isotope'"
                      title="Isotope"
                      regex="int"/>
        <InputElement :input="data.beam.energy"
                      :id="scope+'-measurement-energy'"
                      title="Energy [MeV]"
                      regex="double"/>
        <InputElement :input="data.beam.energy_distribution"
                      :id="scope+'-measurement-energy_distribution'"
                      title="Energy dist [keV]"
                      regex="double" />
        <InputElement :input="data.beam.charge_state"
                      :id="scope+'-measurement-charge_state'"
                      title="Charge state"
                      regex="int" />
      </fieldset>
      <fieldset>
        <legend>Run</legend>
        <InputElement :input="data.beam.spot_size[0]"
                      :id="scope+'-measurement-0'"
                      title="Spot size [mm] X"
                      regex="double" />
        <InputElement :input="data.beam.spot_size[1]"
                      :id="scope+'-measurement-1'"
                      title="Spot size [mm] Y"
                      regex="double" />
<!--        <div class="settings-row">-->
<!--          <span>Spot size [mm]:</span>-->

<!--          <label class="settings-label">x:-->
<!--            <input :id="measurement-0"-->
<!--                   class="settings-input"-->
<!--                   :value="data.beam.spot_size[0]">-->
<!--          </label>-->
<!--          <label class="settings-label">y:-->
<!--            <input :id="measurement-1"-->
<!--                   class="settings-input"-->
<!--                   :value="data.beam.spot_size[1]">-->
<!--          </label>-->
<!--        </div>-->
        <InputElement :input="data.beam.divergence"
                      :id="scope+'-measurement-divergence'"
                      title="Divergence"
                      regex="double" />
        <InputElement :options="['Uniform', 'Gaussian']"
                      :selected="data.beam.profile"
                      :id="scope+'-measurement-profile'"
                      title="Profile"/>
        <InputElement :input="data.run.fluence"
                      :id="scope+'-measurement-fluence'"
                      title="Fluence"
                      regex="double" />
        <InputElement :input="data.run.current"
                      :id="scope+'-measurement-current'"
                      title="Current [nA]"
                      regex="double" />
        <InputElement :input="data.run.charge"
                      :id="scope+'-measurement-charge'"
                      title="Charge [μC]"
                      regex="double" />
        <InputElement :input="data.run.time"
                      :id="scope+'-measurement-time'"
                      title="Time [s]"
                      regex="double" />
      </fieldset>
    </div>
  </b-tab>
</template>

<script>
  import InputElement from '../Elements/InputElement.vue';

  export default {
    name: 'Measurement',
    components: {
      InputElement,
    },
    props: [
      'data',
      'scope',
    ],
    created() { // initialize renamed elements if old settings are used
      this.data.beam.isotope = this.data.beam.ion.replace(/[a-zA-Z]/g, '');
      this.data.beam.element = this.data.beam.ion.replace(/[0-9]/g, '');
    },
  };
</script>

<style>
  #Request-measurement, #Measurement-measurement{
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    flex-wrap: wrap;
  }

  #setup-angle {
    width: 100%;
  }
</style>
