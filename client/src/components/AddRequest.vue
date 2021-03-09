<template>
  <div>
    <h4 v-if="showMessage"  v-b-tooltip.hover
                    title="Yes, it is created!"> {{ message }}</h4>
        
    <b-form @submit="onSubmit" @reset="onReset" class="w-100">
      <b-form-group id="project-name-group">
        <b-form>
          <label for="project-name">Request name</label>
          <b-input v-model="requestName" required id="project-name" type="text"></b-input>
        </b-form>
      </b-form-group>
      <b-button type="submit" variant="primary">Add request</b-button>
      <b-button type="reset" variant="danger">Reset</b-button>
    </b-form>
  </div>
</template>

<script>
  import axios from 'axios';

  export default {
    data() {
      return {
        message: '',
        showMessage: false,
        requestName: '',
      };
    },
    methods: {
      addRequest(payload) {
        axios.post(`http://localhost:5000/request/${payload}`)
          .then(() => {
            this.message = `Request ${payload} created!`;
            this.showMessage = true;
          })
          .catch((error) => {
            // eslint-disable-next-line
            console.error(error);
          });
      },
      // initialize form with blanks
      initForm() {
        this.requestName = '';
      },
      onSubmit(evt) {
        evt.preventDefault();
        const payload = this.requestName;
        this.addRequest(payload);
        this.initForm();
      },
      onReset(evt) {
        evt.preventDefault();
        this.initForm();
      },
    },
    created() {},
  };

</script>
