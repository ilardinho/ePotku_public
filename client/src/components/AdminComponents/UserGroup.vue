<template>
  <!-- Component for showing, editing, creating and deleting user groups -->
  <!-- TODO: Implement functionality -->
  <div class="container">
    <div class="row" >
      <div float>
        <table class="table table-hover">
          <thead>
            <tr>
              <th scope="col">User Group</th>
              <th scope="col">Group ID</th>
              <th scope="col">Group Members</th>
              <th scope="col">Group Info</th>
              <th><button type="button" class="btn btn-success btn-sm">Add User Group</button></th>
            </tr>
          </thead>
          <tbody>
            <tr v-for="(group, index) in groups" :key="index">
              <td>{{ group.userGroupName }}</td>
              <td>{{ group.userGroupID }}</td>
              <td>{{ group.userGroupMembers }}</td>
              <td>{{ group.userGroupInfo }}</td>
              <td>
                <div class="btn-group" role="group">
                  <button type="button" class="btn btn-warning btn-sm">Edit</button>
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

  export default {
    data() {
      return {
        groups: [],
      };
    },
    methods: {
      getUserGroups() {
        const path = 'http://localhost:5000/groups';
        axios.get(path)
          .then((res) => {
            this.groups = res.data.groupdata;
          })
          .catch((error) => {
            // eslint-disable-next-line
            console.error(error);
          });
      },
    },
    created() {
      this.getUserGroups();
    },
  };
</script>
