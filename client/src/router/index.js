import Vue from 'vue';
import VueRouter from 'vue-router';
import FrontPage from '../views/FrontPage.vue';
import Application from '../views/Application.vue';
import Projects from '../components/Projects.vue';
import Admin from '../views/AdminConsole.vue';

Vue.use(VueRouter);

const routes = [
  {
    path: '/',
    name: 'FrontPage',
    component: FrontPage,
  },
  {
    path: '/app',
    name: 'App',
    component: Application,
  },
  {
    path: '/admin',
    name: 'Admin',
    component: Admin,
  },
  {
    path: '/projects',
    name: 'Projects',
    component: Projects,
  },
];

const router = new VueRouter({
  mode: 'history',
  base: process.env.BASE_URL,
  routes,
});

export default router;
