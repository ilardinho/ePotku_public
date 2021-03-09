<template>
  <b-modal id="newMeasurement"  hide-footer title="Add a new measurement">
     <h4 v-if="showMessage"  v-b-tooltip.hover
                    title="Yes, this happened!"> {{ message }}</h4>
    <!-- TODO: Remove this. Request where measurement will be uploaded
    should be get from VUEX, not by typing it here. -->
    <b-form-group id="request-name-group">
      <b-form>
        <label for="request-name">Request name</label>
        <b-input v-model="request_name" id="request-name" required></b-input>
      </b-form>
    </b-form-group>
    <!-- Name field for sample -->
    <b-form-group id="sample-name-group">
      <b-form>
        <label for="sample-name">Sample name</label>
        <b-input v-model="sample_name" id="sample-name" required></b-input>
      </b-form>
    </b-form-group>
    <!-- Name field for measurement -->
    <b-form-group id="measurement-name-group">
      <b-form>
        <label for="measurement-name">Measurement name</label>
        <b-input v-model="measurement_name" id="measurement-name" required></b-input>
      </b-form>
    </b-form-group>
    <!-- Input field for measurement file-->
    <b-form-group>
      <b-form-file v-model="file" ref="file"  :state="Boolean(file)"
        v-on:change="onChangeFileUpload()" placeholder="Choose a file or drop it here..."
        drop-placeholder="Drop file here..." enctype="multipart/form-data">
      </b-form-file>
      <div class="mt-3">Selected file: {{ file ? file.name : '' }}</div>
    </b-form-group>
    <!-- <b-button type="submit" variant="primary" v-on:click="submitFile(file)">
      Upload</b-button> -->
    <b-button variant="primary" v-on:click="submitFile(file)">Upload</b-button>
  </b-modal>
</template>

<script>
  import axios from 'axios';

  export default {
    name: 'NewMeasurement',
    data() {
      return {
        message: '',
        showMessage: false,
        file: null,
        sample_name: '',
        measurement_name: '',
        request_name: '',
      };
    },
    methods: {
      submitFile() {
        const formData = new FormData();
        const payload = this.request_name;
        formData.append('file', this.file);
        formData.append('sample_name', this.sample_name);
        formData.append('measurement_name', this.measurement_name);
        formData.append('request_name', this.request_name);
        axios.post(`http://localhost:5000/measurement/upload/${payload}`, formData, {
          headers: {
            'Content-Type': 'multipart/form-data',
          },
        })
          .then(() => {
            this.message = `Measurement ${this.measurement_name} added!`;
            this.showMessage = true;
            this.initForm();
          })
          .catch((error) => {
            this.message = `Something went wrong when adding ${this.measurement_name}! ${this.error}`;
            this.showMessage = true;
            // eslint-disable-next-line
            console.error(error);
          });
      },
      // initialize form with blanks
      initForm() {
        this.file = null;
        this.sample_name = '';
        this.measurement_name = '';
        this.request_name = '';
      },
      // onSubmit(evt) {
      //   evt.preventDefault();
      //   this.submitFile();
      //   this.$refs[this.newMeasurement].hide();
      // },
      onChangeFileUpload() {
        // const j = 0;
        // this.file = this.$refs.file.files[j];
        // const j = 0;
        this.file = this.$refs.file.file;
      },
    },
    created() {},
  };
</script>
