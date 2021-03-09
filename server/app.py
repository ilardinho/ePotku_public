# coding=utf-8
"""
ePotku is a web application version of Potku, 
a graphical user interface for analyzation and
visualization of measurement data collected from a ToF-ERD
telescope. For physics calculations Potku uses external
analyzation components.
Copyright (C) 2020 Minja Hänninen, Ilari Jalli, Ville Kuokkanen,
Pasi Niininen and Tuomas Pitkänen
​
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
​
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
​
You should have received a copy of the GNU General Public License
along with this program (file named 'LICENCE').
"""
__author__ = "Minja Hänninen \n Ilari Jalli \n Ville Kuokkanen \n" \
             "Pasi Niininen \n ""Tuomas Pitkänen"
__version__ = "1.0"

# -*- coding: utf-8 -*-
import functools
# import logging
import os
import json
import sys
import time

from flask import Flask, jsonify, request, json, flash, redirect
# from flask_caching import Cache
from flask_cors import CORS
from flask_jwt_extended import (
    JWTManager, jwt_required, get_jwt_identity,
    create_access_token, verify_jwt_in_request
)

import call_interface
from models import db, User

# ALLOWED_EXTENSIONS = {'txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif', 'asc'}
ALLOWED_EXTENSIONS = {'asc'}

# TODO: Replace status messages with status codes

# Cache-related stuff:
# config = {
#     "CACHE_TYPE": "simple", # Flask-Caching related configs
#     "CACHE_DEFAULT_TIMEOUT": 3600
# }

# instantiate the app
app = Flask(__name__)

# Settings
app.config.from_pyfile('default_settings.cfg')
# TODO: Skip envvar if it is not set (silent=True)?
#       Display a warning in that case.
#       Do not warn in development mode
app.config.from_envvar('EPOTKU_SETTINGS', silent=True)

# cache = Cache(app)

# enable CORS
CORS(app, resources={r'/*': {'origins': '*'}})

# Avoid circular dependencies by initializing SQLAlchemy with db from models.py
db.init_app(app)
# TODO: Check for existing tables before recreating
with app.app_context():
    db.create_all()

# Initialize Flask-JWT-Extended
jwt = JWTManager(app)

# TODO: Remove test users
# Add test users
with app.app_context():
    if User.query.filter(User.email == 'user@example.com').first() is None:
        test_user = User(email='user@example.com', first_name='User',
                         last_name='Test')
        test_user.set_password('Password1')
        db.session.add(test_user)
        db.session.commit()
    if User.query.filter(User.email == 'admin@example.com').first() is None:
        test_admin = User(email='admin@example.com', first_name='Admin',
                          last_name='Istrator', is_admin=True)
        test_admin.set_password('Password1')
        db.session.add(test_admin)
        db.session.commit()


# TODO: Move to another file
def admin_required(func):
    """Wrapper for denying access to an endpoint by non-admin users."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # TODO: Check claims in JWT instead (without database):
        # https://flask-jwt-extended.readthedocs.io/en/stable/custom_decorators/
        # TODO: Check for JWT freshness (verify_fresh_jwt_in_request)
        verify_jwt_in_request()
        current_user_id = get_jwt_identity()
        user = User.query.get(current_user_id)
        if not user.is_admin:
            return jsonify({'message': 'Admin role required'}), 403
        else:
            return func(*args, **kwargs)

    return wrapper


# TODO: get settings from the appropriate request folder
@app.route('/settings/Request', methods=['GET'])
def get_request_settings():
    settings = {
        'measurement': None,
        'detector': None,
        'profile': None,
        'mcsimu': None
    }
    for setting in settings:
        fp = os.path.join(app.static_folder, 'request_settings',
                          'Default.' + setting)
        with open(fp) as file:
            settings[setting] = json.load(file)
    return json.dumps(settings)


# TODO: get settings from the appropriate request folder
@app.route('/settings/Measurement', methods=['GET'])
def get_measurement_settings():
    settings = {
        'measurement': None,
        'detector': None,
        'profile': None,
    }
    for setting in settings:
        fp = os.path.join(app.static_folder, 'measurement_settings',
                          'Default.' + setting)
        with open(fp) as file:
            settings[setting] = json.load(file)
    return settings


@app.route('/login', methods=['POST'])
def login():
    """Handle user login. Return an access token if login was successful."""
    if not request.is_json:
        return jsonify({'message': 'Missing JSON in request'}), 400

    # Check for missing parameters
    email = request.json.get('email')
    password = request.json.get('password')
    if email is None:
        return jsonify({'message': 'Missing email parameter'}), 400
    if password is None:
        return jsonify({'message': 'Missing password parameter'}), 400

    # TODO: Token expiry?
    # Check credentials
    # Assumption: email is unique (email is set as unique in database)
    user = User.query.filter(User.email == email).first()
    if user is not None and user.check_password(password):
        access_token = create_access_token(identity=user.id)
        user_dict = user.get_dict_without_password()
        return jsonify({'access_token': access_token, 'user': user_dict}), 200
    return jsonify({'message': 'Bad username or password'}), 401


# TODO: Remove test page
@app.route('/jwt_protected')
@jwt_required
def jwt_protected():
    current_user_id = get_jwt_identity()
    user = User.query.get(current_user_id)
    return jsonify({'user': user.id}), 200


# TODO: Remove test page
@admin_required
def admin_protected():
    current_user_id = get_jwt_identity()
    user = User.query.get(current_user_id)
    return jsonify({'user': user.id, 'admin': user.is_admin}), 200


# TODO: Use flask.views for REST methods?


# TODO: Add pagination?
@app.route('/user', methods=['GET'])
@admin_required
def get_users():
    """Handler for getting all users' information.

    Does not return passwords.
    """
    users_dict = [user.get_dict_without_password()
                  for user in User.query.all()]
    return jsonify({'users': users_dict}), 200


# TODO: Need to handle expiration. Something like {if no input from user for x minutes: clear()}
# TODO: Use cache for big files
# def set_cache():
#     cache.set('caching_data', True) # for requests received while caching in progress
#     # cached data is a 2d array from .asc in potku/modules/measurement.py load_data()
#     data = call_interface.get_request_measurement()
#     cache.set('data', data)
#     cache.delete('caching_data')
#     return 'caching successful'
#
#
# def check_cache():
#     if cache.get('data') is None: # data not in cache
#         if cache.get('caching_data') is None: # caching not in progress
#             set_cache()
#         else:
#             n=0
#             # TODO: more sophisticated timeout handling
#             while cache.get('data') is None and n < 30: # wait for caching to finish
#                 time.sleep(1)
#                 n+=1
#     return cache.get('data')


@app.route('/tofe', methods=['GET'])
def get_tofe_data():
    # data = check_cache()
    compression_x = int(request.args['compressionX'])
    compression_y = int(request.args['compressionY'])
    return call_interface.load_testrequest(compression_x, compression_y)


@app.route('/selections', methods=['GET'])
def get_selections():
    # data = check_cache()
    return jsonify(call_interface.load_selections())


@app.route('/user/<int:user_id>', methods=['GET'])
@admin_required
def get_user(user_id):
    """Handler for getting a user's information.

    Does not return passwords.
    """
    # TODO: Add error handling
    user = User.query.get(user_id)
    if user is None:
        return '', 404
    return jsonify({'user': user.get_dict_without_password()}), 200


@app.route('/user', methods=['POST'])
@admin_required
def add_user():
    """Handler for creating a new user"""
    # TODO: Add error handling
    # FIXME: Currently returns HTTP 500 if email is not unique.
    #        Return 409 Conflict instead?
    # TODO: Parse user information in a separate function
    email = request.json.get('email')
    password = request.json.get('password')
    first_name = request.json.get('first_name')
    last_name = request.json.get('last_name')
    is_admin = request.json.get('is_admin')
    new_user = User(email=email, first_name=first_name, last_name=last_name,
                    is_admin=is_admin)
    new_user.set_password(password)
    db.session.add(new_user)
    db.session.commit()
    return jsonify(new_user.get_dict_without_password()), 201


# A PUT method for users is not defined because:
# - passwords should not be modified this way (use another endpoint for
#   that). Administrators do not know their users' passwords.
# - database handles setting the id for new users, not the end user.
# - PATCH fits better.


@app.route('/user/<int:user_id>', methods=['PATCH'])
@admin_required
def update_user(user_id):
    """Handler for updating a user's information.

    Passwords cannot be changed this way.
    """
    # TODO: Add error handling
    updated_user = User.query.get(user_id)
    email = request.json.get('email')
    first_name = request.json.get('first_name')
    last_name = request.json.get('last_name')
    is_admin = request.json.get('is_admin')
    # (Do not change password)
    if email is not None:
        updated_user.email = email
    if first_name is not None:
        updated_user.first_name = first_name
    if last_name is not None:
        updated_user.last_name = last_name
    if is_admin is not None:
        updated_user.is_admin = is_admin
    db.session.commit()
    return jsonify(updated_user.get_dict_without_password()), 200


@app.route('/user/<int:user_id>', methods=['DELETE'])
@admin_required
def delete_user(user_id):
    """Handler for deleting a user."""
    if user_id is None:
        return jsonify({'message': 'Missing user_id query parameter'}), 404

    user = User.query.get(user_id)
    if user is None:
        # DELETE is idempotent, therefore allow "deleting" a nonexistent user.
        return '', 204
    db.session.delete(user)
    db.session.commit()
    return '', 204


# TODO: logging for the application
# Logging test, requires 'import logging'. This writes logs
# which is required info for admins
# logging.basicConfig(filename='loki.json',
# level=logging.DEBUG,
# format='%(asctime)s %(levelname)s %(name)s %(threadName)s : %(message)s')


@app.route('/cuts', methods=['GET'])
def get_cuts():
    # data = check_cache()
    return call_interface.get_cuts()


@app.route('/cuts', methods=['POST'])
def save_cuts():
    # data = check_cache()
    return call_interface.save_cuts((request.get_json())['cuts'])


@app.route('/depth_profiles', methods=['GET'])
def get_depth_profiles():
    # data = check_cache()
    depth_units = request.args['xUnits']
    cuts = request.args.getlist('selectedCuts[]')
    error = float(request.args['systematicError'])
    ref_density = float(request.args['referenceDensity'])
    return call_interface.load_depth(depth_units, cuts, error)


@app.route('/element_losses', methods=['GET'])
def get_element_losses():
    # data = check_cache()
    splits = int(request.args['splits'])
    ref_cut = request.args['referenceCut']
    cuts = request.args.getlist('cuts[]')
    y_axis = request.args['yAxis']
    return call_interface.load_element_losses(cuts, ref_cut, splits)


@app.route('/element_losses', methods=['POST'])
def save_element_losses():
    # data = check_cache()
    args = (request.get_json())
    splits = int(args['splits'])
    ref_cut = args['referenceCut']
    cuts = args['cuts']
    return call_interface.save_splits(cuts, ref_cut, splits)


@app.route('/energy_spectra', methods=['GET'])
def get_energy_spectra():
    # data = check_cache()
    cuts = request.args.getlist('cuts[]')
    bin_width = float(request.args['binWidth'])
    return call_interface.load_energy_spectra(cuts, bin_width)


@app.route('/request/<requestName>', methods=['POST'])
def add_request(requestName):
    """Handler for creating a new request"""
    if requestName is None:
        return jsonify({'message': 'Missing requestName query parameter'}), 404
    response_object = {'status': 'success'}
    # TODO: request owner UUID should be provided within HTTP request
    requestOwner = '338d36aa-9f96-11ea-b8a4-b8e8560bf606'
    published = '1'
    response_object['request'] = call_interface.new_request(requestName,
                                                            requestOwner,
                                                            published)
    return jsonify(response_object), 201


# TODO: change endpoint to '/request'
@app.route('/request/<requestName>', methods=['DELETE'])
def delete_request(requestName):
    """Handler for deleting a request"""
    if requestName is None:
        return jsonify({'message': 'Missing requestName query parameter'}), 404
    response_object = {'status': 'success'}
    call_interface.delete_request(requestName)
    response_object['Info'] = 'Request deleted'
    return jsonify(response_object), 200


# TODO: change endpoint to '/request?user_id=xyz'
@app.route('/allrequests/<ownerUUID>', methods=['GET'])
def get_requests(ownerUUID):
    """Handler for getting users requests"""
    if ownerUUID is None:
        return jsonify({'message': 'Missing ownerUUID query parameter'}), 404
    allrequests = call_interface.get_requests(ownerUUID)
    response_object = {'status': 'success', 'requests': allrequests}
    return response_object, 200


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


# TODO: Change endpoint to '/request/<request_id>'
@app.route('/requestcontent/<request_name>', methods=['GET'])
def get_request(request_name):
    """Handler for getting request contents"""
    if request_name is None:
        return jsonify({'message': 'Missing request_name query parameter'}), 404
    request_content = call_interface.request_content(request_name)
    response_object = {'status': 'success', 'samples': request_content}
    print(response_object)
    return response_object, 200


# TODO: This was called upload_test. Is something unfinished?
# TODO: Change endpoint to '/request/<request_id>/measurement'
@app.route('/measurement/upload/<requestName>', methods=['POST'])
def add_measurement(requestName):
    """Handler for uploading a new measurement"""
    sample_name = request.form['sample_name']
    measurement_name = request.form['measurement_name']
    if 'file' not in request.files:
        flash('No file part')
        return "There's no file to upload", 400
    measurement_file = request.files['file']
    if measurement_file.filename == '':
        flash('No selected file')
        return "There's no file to upload", 400
    if measurement_file and allowed_file(measurement_file.filename):
        upload_measurement = call_interface.measurement_upload(measurement_file,
                                                              requestName,
                                                              sample_name,
                                                              measurement_name)
        response_object = {'status': 'success',
                           'Created': upload_measurement}
        print(response_object)
        return response_object, 201
    return 'File is not supported, please upload .asc file!', 415


if __name__ == '__main__':
    app.run()
