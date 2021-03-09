/* eslint-disable no-param-reassign */
/* Rule conflicts with Vuex. */

import Vue from 'vue';
import Vuex from 'vuex';
import VuexPersistence from 'vuex-persist';

import axios from 'axios';

const BACKEND_URL = 'http://localhost:5000/';
const LOGIN_URL = `${BACKEND_URL}login`;

Vue.use(Vuex);

// TODO: Move authentication stuff to its own service and then use it here
export default new Vuex.Store({
  state: {
    user: {
      accessToken: '',
      // TODO: Give a better name:
      attributes: {},
    },
  },
  plugins: [
    new VuexPersistence().plugin,
  ],
  getters: {
    getAuthHeader(state, getters) {
      // TODO: Throw exception?
      if (!getters.isAuthenticated) return '';
      return `Bearer ${state.user.accessToken}`;
    },
    isAuthenticated(state) {
      // Convert result to boolean
      return !!(state.user.accessToken && state.user.accessToken.length > 0);
    },
  },
  actions: {
    login(context, payload) {
      return axios.post(LOGIN_URL, payload.credentials)
        .then((response) => {
          context.commit('setAccessToken', response.data.access_token);
          context.commit('setUser', response.data.user);
          // eslint-disable-next-line no-console
          console.log(response);
        })
        .catch((error) => {
          // eslint-disable-next-line no-console
          console.error(error);
        });
    },
    logout(context) {
      context.commit('clearAccessToken');
      context.commit('clearUser');
    },
  },
  mutations: {
    setAccessToken(state, token) {
      state.user.accessToken = token;
    },
    clearAccessToken(state) {
      state.user.accessToken = '';
    },
    setUser(state, attributes) {
      state.user.attributes = attributes;
    },
    clearUser(state) {
      state.user.attributes = {};
    },
  },
});
