<template>
  <!-- This container is to show the projects. This one makes GET query to flask-server -->
  <div class="container">
    <div class="row">
      <div>
        <h4 v-if="showMessage" v-b-tooltip.hover title="Yes, this happened"> {{ message }}</h4>

        <table class="table table-hover">
          <thead>
            <tr>
              <th scope="col">Request name</th>
              <th scope="col">Request ID</th>
              <th scope="col">Published</th>
              <th></th>
            </tr>
          </thead>
          <tbody>
            <tr v-for="(request, index) in requests" :key="index">
              <td>{{ request.requestName }}</td>
              <td>{{ request.requestUUID }}</td>
              <td>
                <b-icon-circle-fill v-if="request.published == 1" variant="success" font-scale="2">
                </b-icon-circle-fill>
              </td>
              <td>
                <div class="btn-group" role="group">
                  <button type="button" class="btn btn-success btn-sm"
                    @click="openRequest(request.requestName)" v-b-tooltip.hover
                    title="Open request">Open</button>

                  <!-- TODO: Publishing is not done yet -->
                  <button type="button" class="btn btn-warning btn-sm"
                    @click="publishRequest(request.requestName)" v-b-tooltip.hover
                    title="Publish request as open data">Publish</button>

                  <button type="button" class="btn btn-danger btn-sm"
                    @click="onDeleteRequest(request)" v-b-tooltip.hover
                    title="Delete request">Delete</button>
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
        message: '',
        showMessage: false,
        requests: [],
        request_content: [],
      };
    },
    methods: {
      getRequests(payload) {
        const path = `http://localhost:5000/allrequests/${payload}`;
        axios.get(path)
          .then((res) => {
            this.requests = res.data.requests;
          })
          .catch((error) => {
            // eslint-disable-next-line
            console.error(error);
          });
      },
      removeRequest(requestName) {
        const path = `http://localhost:5000/request/${requestName}`;
        axios.delete(path)
          .then(() => {
            // TODO: getRequest should be called with proper userUUID
            this.getRequests('338d36aa-9f96-11ea-b8a4-b8e8560bf606');
            this.message = `Request ${requestName} deleted!`;
            this.showMessage = true;
          })
          .catch((error) => {
            // eslint-disable-next-line
            console.error(error);
            this.getRequests();
          });
      },
      openRequest(requestName) {
        const path = `http://localhost:5000/requestcontent/${requestName}.potku`;
        axios.get(path)
          .then((res) => {
            // TODO: getRequest should be called with proper userUUID
            this.getRequests('338d36aa-9f96-11ea-b8a4-b8e8560bf606');
            this.message = `Request ${requestName} opened!`;
            this.request_content = res.data.requests;
            this.showMessage = true;
          })
          .catch((error) => {
            // eslint-disable-next-line
            console.error(error);
            this.getRequests();
          });
      },
      onDeleteRequest(request) {
        // const r = confirm(`Are you sure you want to delete ${request} request? `);
        if (window.confirm(`Are you sure you want to delete ${request.requestName} request? `)) {
          this.removeRequest(request.requestName);
        }
      },
    },
    created() {
      const payload = '338d36aa-9f96-11ea-b8a4-b8e8560bf606';
      this.getRequests(payload);
    },
  };

</script>
