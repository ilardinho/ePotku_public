<template>
  <!-- Modal for deleting a user -->
  <b-modal id="userDeleteModal" :title="'Delete user: ' + user.email + '?'"
    ok-title="Delete" ok-variant="danger"
    @ok="onSubmit" @cancel="onCancel" @hide="onHide"
  >
    <UserForm id="userForm" :user="this.user" :isDisabled="true" />
  </b-modal>
</template>

<script>
  import axios from 'axios';

  import UserForm from './UserForm.vue';

  const BACKEND_URL = 'http://localhost:5000/';
  const USER_URL = `${BACKEND_URL}user/`;

  export default {
    name: 'UserDeleteModal',
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
        this.deleteUserData();
        this.clearForm();
      },
      onCancel() {
        this.clearForm();
      },
      onHide() {
        this.clearForm();
      },

      /* HTTP methods */

      // TODO: Reduce copy-paste from Users, UserCreateModal, UserEditModal and UserDeleteModal
      //       with a new UserService javascript file that handles different HTTP methods.

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
      deleteUserData() {
        axios.delete(this.getCurrentUserUrl(), {
          headers: {
            Authorization: this.$store.getters.getAuthHeader,
          },
        })
          .then((response) => {
            // eslint-disable-next-line no-console
            console.log(response);
            // Notify parent about user deletion
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
