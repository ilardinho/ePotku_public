<template>
  <b-modal id="globalSettings" title="Global settings" ok-title="Save"
           no-close-on-backdrop no-close-on-esc size="lg" scrollable>
    <b-tabs>
      <b-tab id="tof-histogram" title="ToF-E Histogram">
        <InputElement :options="['Default color','Greyscale','Greyscale (inverted)']"
                      selected="Default color"
                      title="Color scheme" />

        <input type="checkbox" id="invert-x"><label for="invert-x">Invert X axis</label><br>
        <input type="checkbox" id="invert-y"><label for="invert-y">Invert Y axis</label><br>
        <input type="checkbox" id="transpose"><label for="transpose">Transpose axes</label><br>


        <InputElement input=""
                      title="Compression [X axis]" />
        <InputElement input=""
                      title="Compression [Y axis]" />
        <InputElement :options="['Automatic', 'Manual']"
                      selected="Automatic"
                      title="Range" />
      </b-tab>
      <b-tab id="depth-profiles" title="Depth Profiles">
        <InputElement input=""
                      title="Number if iterations"
                      />
        <InputElement :options="['Rutherford', 'L\'ecuyer', 'Andersen']"
                      selected="Rutherford"
                      title="Cross section"/>
      </b-tab>
      <b-tab id="import" title="Import">
        <span class="settings-row">ADC 0
          <InputElement input=""
                        title="Low" />
          <InputElement input=""
                        title="High" />
        </span>
          <span class="settings-row">ADC 1
          <InputElement input=""
                        title="Low" />
          <InputElement input=""
                        title="High" />
        </span>
        <span class="settings-row">ADC 2
          <InputElement input=""
                        title="Low" />
          <InputElement input=""
                        title="High" />
        </span>
        <InputElement input=""
                      title="Preview count" />
      </b-tab>
    </b-tabs>
  </b-modal>
</template>


<script>
  import InputElement from './Elements/InputElement.vue';

  export default {
    components: {
      InputElement,
    },
  };
</script>

<style>

</style>
