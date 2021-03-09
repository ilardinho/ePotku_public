<style scoped>
  .sidebar a {
    padding: 8px 8px 8px 32px;
    text-decoration: none;
    color: #000000;
    display: block;
    transition: 0.3s;
  }

  .sidebar a:hover {
    color: #515050;
  }

</style>

<template>
  <div class="no-overflow">
    <!--Left navbar -->
    <div class="bg-light sidebar">
      <!-- This is the navbar part -->
      <b-navbar-nav>
        <ul class="nav flex-column">

          <!-- navigation -->
          <b-nav vertical align="left">
            <b-nav-item v-b-modal.Request-settings>
              <b-icon icon="gear" font-scale="1.5"></b-icon> Request Settings
            </b-nav-item>
            <b-nav-item v-b-modal.Measurement-settings>
              <b-icon icon="gear-fill" font-scale="1.5"></b-icon> Measurement Settings
            </b-nav-item>
          </b-nav>

          <!--Project card -->
          <!-- TODO: max width for navbar, wrapping for text -->
          <b-card header="Request: Demo request">
            <b-list-group>
              <b-list-group-item>
                <b-link @click="openRequest('Potku_v2')">
                  <b-icon icon="folder" font-scale="1.5"></b-icon> Potku_v2
                </b-link>
              </b-list-group-item>

              <b-list-group-item>
                <table class="table table-hover">
                  <tbody v-for="(sample, index) in request_content" :key="index">

                    <tr class="name-row">
                      <th class="name" colspan="3">{{index}}</th>
                    </tr>
                    <tr v-for="(sample, index) in sample" :key="index">
                      <td class="joo" colspan="3">{{sample}}</td>
                    </tr>
                  </tbody>
                </table>
              </b-list-group-item>
            </b-list-group>
          </b-card>

        </ul>
      </b-navbar-nav>
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
      openRequest(requestName) {
        const path = `http://localhost:5000/requestcontent/${requestName}.potku`;
        axios.get(path)
          .then((res) => {
            // TODO: getRequest should be called with proper userUUID
            this.getRequests('338d36aa-9f96-11ea-b8a4-b8e8560bf606');
            this.message = `Request ${requestName} opened!`;
            this.request_content = res.data.samples;
            this.showMessage = true;
          })
          .catch((error) => {
            // eslint-disable-next-line
            console.error(error);
            this.getRequests();
          });
      },
      onDeleteRequest(request) {
        this.removeRequest(request.requestName);
      },
    },
    created() {
      const payload = '338d36aa-9f96-11ea-b8a4-b8e8560bf606';
      this.getRequests(payload);
      this.openRequest('Potku_v2');
    },
  };

</script>
