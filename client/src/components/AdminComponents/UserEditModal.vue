<template>
  <!-- Modal for editing an exiting user's information -->
  <b-modal id="userEditModal" :title="'Editing user: ' + user.email"
    ok-title="Save changes"
    @ok="onSubmit" @cancel="onCancel" @hide="onHide"
  >
    <UserForm id="userForm" :user="this.user" :isDisabled="false" />
  </b-modal>
</template>

<script>
  import axios from 'axios';

  import UserForm from './UserForm.vue';

  const BACKEND_URL = 'http://localhost:5000/';
  const USER_URL = `${BACKEND_URL}user/`;

  export default {
    name: 'UserEditModal',
    components: {
      UserForm,
    },
    data() {
      return {
        user: {},
        selectedUserId: null,
      };
    },
    methods: {
      loadForm(userId) {
        this.selectedUserId = userId;
        this.getUserData();
      },
      clearForm() {
        this.user = {};
        this.selectedUserId = null;
      },
      onSubmit() {
        this.patchUserData();
        this.clearForm();
      },
      onCancel() {
        this.clearForm();
      },
      onHide() {
        this.clearForm();
      },

      /* HTTP methods */

      getCurrentUserUrl() {
        return `${USER_URL}${this.selectedUserId}`;
      },

      // TODO: Display "success" and "failure" messages for HTTP methods

      getUserData() {
        axios.get(this.getCurrentUserUrl(), {
          headers: {
            Authorization: this.$store.getters.getAuthHeader,
          },
        })
          .then((response) => {
            this.user = response.data.user;
          })
          .catch((error) => {
            // eslint-disable-next-line no-console
            console.error(error);
          });
      },
      patchUserData() {
        axios.patch(this.getCurrentUserUrl(), this.user, {
          headers: {
            Authorization: this.$store.getters.getAuthHeader,
          },
        })
          .then((response) => {
            // eslint-disable-next-line no-console
            console.log(response);
            // Notify parent about user update
            this.$emit('update');
          })
          .catch((error) => {
            // eslint-disable-next-line no-console
            console.error(error);
          });
      },
    },
  };
</script>
