<template>
  <!-- Component for showing, editing, creating and deleting users -->
  <div class="container">
    <UserCreateModal ref="userCreateModal" @update="getUsers" />
    <UserEditModal ref="userEditModal" @update="getUsers" />
    <UserDeleteModal ref="userDeleteModal" @update="getUsers" />
    <div class="row" >
      <div float>
        <table class="table table-hover">
          <thead>
            <tr>
              <!-- TODO: Is logged in status feasible to implement (especially with JWT)? -->
              <!-- <th scope="col">Logged in</th> -->
              <th scope="col">Name</th>
              <th scope="col">Username (TODO)</th>
              <th scope="col">Email</th>
              <th scope="col">User groups (TODO)</th>
              <th scope="col">Is admin</th>
              <th scope="col">Priority (TODO)</th>
              <th scope="col">User ID</th>
              <th>
                <b-button
                  v-b-modal.userCreateModal
                  class="btn btn-success btn-sm"
                  @click="createUser()"
                >
                  Create a user
                </b-button>
              </th>
            </tr>
          </thead>
          <tbody>
            <tr v-for="(user, index) in users" :key="index">
              <!-- <td><b-icon-circle-fill variant="success"></b-icon-circle-fill></td> -->
              <td>{{ user.first_name + ' ' + user.last_name}}</td>
              <td>{{ user.username }}</td>
              <td>{{ user.email }}</td>
              <td>{{ user.userGroup }}</td>
              <td>{{ user.is_admin }}</td>
              <td>{{ user.priority }}</td>
              <td>{{ user.id }}</td>
              <td>
                <div class="btn-group" role="group">
                  <b-button
                    v-b-modal.userEditModal
                    class="btn btn-warning btn-sm"
                    @click="editUser(user.id)"
                  >
                    Edit
                  </b-button>
                  <b-button
                    v-b-modal.userDeleteModal
                    class="btn btn-danger btn-sm"
                    @click="deleteUser(user.id)"
                  >
                    Delete
                  </b-button>
                </div>
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>
  </div>
</template>

<script>
  import axios from 'axios';

  import UserCreateModal from './UserCreateModal.vue';
  import UserEditModal from './UserEditModal.vue';
  import UserDeleteModal from './UserDeleteModal.vue';

  const BACKEND_URL = 'http://localhost:5000';
  const USER_URL = `${BACKEND_URL}/user`;

  export default {
    name: 'User',
    components: {
      UserCreateModal,
      UserEditModal,
      UserDeleteModal,
    },
    data() {
      return {
        users: [],
      };
    },
    methods: {
      createUser() {
        this.$refs.userCreateModal.loadForm();
      },
      editUser(userId) {
        this.$refs.userEditModal.loadForm(userId);
      },
      deleteUser(userId) {
        this.$refs.userDeleteModal.loadForm(userId);
      },

      /* HTTP methods */

      getUsers() {
        axios.get(USER_URL, {
          headers: {
            Authorization: this.$store.getters.getAuthHeader,
          },
        })
          .then((res) => {
            this.users = res.data.users;
          })
          .catch((error) => {
            // TODO: Display an error message
            // eslint-disable-next-line no-console
            console.error(error);
          });
      },
    },
    created() {
      this.getUsers();
    },
  };
</script>
