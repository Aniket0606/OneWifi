/************************************************************************************
  If not stated otherwise in this file or this component's Licenses.txt file the
  following copyright and licenses apply:

  Copyright 2018 RDK Management

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <pthread.h>
#include <unistd.h>
#include "collection.h"
#include "wifi_webconfig.h"
#include "wifi_monitor.h"
#include "wifi_util.h"
#include "wifi_ctrl.h"

#ifdef EM_APP
webconfig_subdoc_object_t   em_channel_stats_objects[3] = {
    { webconfig_subdoc_object_type_version, "Version" },
    { webconfig_subdoc_object_type_subdoc, "SubDocName" },
    { webconfig_subdoc_object_type_stats_config, "ChannelScanResponse" },
};

webconfig_error_t init_em_channel_stats_subdoc(webconfig_subdoc_t *doc)
{
    doc->num_objects = sizeof(em_channel_stats_objects)/sizeof(webconfig_subdoc_object_t);
    memcpy((unsigned char *)doc->objects, (unsigned char *)&em_channel_stats_objects, sizeof(em_channel_stats_objects));

    return webconfig_error_none;
}

webconfig_error_t access_check_em_channel_stats_subdoc(webconfig_t *config, webconfig_subdoc_data_t *data)
{
    return webconfig_error_none;
}

webconfig_error_t translate_from_em_channel_stats_subdoc(webconfig_t *config, webconfig_subdoc_data_t *data)
{
    if (((data->descriptor & webconfig_data_descriptor_translate_to_ovsdb) == webconfig_data_descriptor_translate_to_ovsdb) ||  ((data->descriptor & webconfig_data_descriptor_translate_to_easymesh) == webconfig_data_descriptor_translate_to_easymesh)) {
        if (config->proto_desc.translate_to(webconfig_subdoc_type_em_channel_stats, data) != webconfig_error_none) {
            if ((data->descriptor & webconfig_data_descriptor_translate_to_ovsdb) == webconfig_data_descriptor_translate_to_ovsdb) {
                return webconfig_error_translate_to_ovsdb;
            } else {
                return webconfig_error_translate_to_easymesh;
            }
        }
    } else if ((data->descriptor & webconfig_data_descriptor_translate_to_tr181) == webconfig_data_descriptor_translate_to_tr181) {

    } else {
        // no translation required
    }
    //no translation required
    return webconfig_error_none;
}

webconfig_error_t translate_to_em_channel_stats_subdoc(webconfig_t *config, webconfig_subdoc_data_t *data)
{
    if (((data->descriptor & webconfig_data_descriptor_translate_from_ovsdb) == webconfig_data_descriptor_translate_from_ovsdb) || ((data->descriptor & webconfig_data_descriptor_translate_from_easymesh) == webconfig_data_descriptor_translate_from_easymesh)) {
        if (config->proto_desc.translate_from(webconfig_subdoc_type_em_channel_stats, data) != webconfig_error_none) {
            if ((data->descriptor & webconfig_data_descriptor_translate_from_ovsdb) == webconfig_data_descriptor_translate_from_ovsdb) {
                return webconfig_error_translate_from_ovsdb;
            } else {
                return webconfig_error_translate_from_easymesh;
            }
        }
    } else if ((data->descriptor & webconfig_data_descriptor_translate_from_tr181) == webconfig_data_descriptor_translate_from_tr181) {

    } else {
        // no translation required
    }
    return webconfig_error_none;
}

webconfig_error_t encode_em_channel_stats_subdoc(webconfig_t *config, webconfig_subdoc_data_t *data)
{
    cJSON *json, *st_obj_arr;
    webconfig_subdoc_decoded_data_t *params;
    char *str;
    mac_addr_str_t mac_str;

    if (data == NULL) {
        wifi_util_error_print(WIFI_WEBCONFIG, "%s:%d: NULL data Pointer\n", __func__, __LINE__);
        return webconfig_error_encode;
    }

    params = &data->u.decoded;
    if (params == NULL) {
        wifi_util_error_print(WIFI_WEBCONFIG, "%s:%d: NULL Pointer\n", __func__, __LINE__);
        return webconfig_error_encode;
    }

    json = cJSON_CreateObject();
    if (json == NULL) {
        wifi_util_error_print(WIFI_WEBCONFIG, "%s:%d: json create object failed\n", __func__,
            __LINE__);
        return webconfig_error_encode;
    }

    data->u.encoded.json = json;

    cJSON_AddStringToObject(json, "Version", "1.0");
    cJSON_AddStringToObject(json, "SubDocName", "EM_Channel_Stats");

    channel_scan_response_t *neigh_stats = params->collect_stats.stats;

    if ((neigh_stats == NULL) || (neigh_stats->num_results == 0)) {
        wifi_util_dbg_print(WIFI_WEBCONFIG, "%s:%d: No neighbor stats to encode \n", __func__,
            __LINE__);
    }

    to_mac_str(neigh_stats->ruid, mac_str);
    cJSON_AddStringToObject(json, "ScannerMac", mac_str);

    st_obj_arr = cJSON_CreateArray();
    cJSON_AddItemToObject(json, "ChannelScanResponse", st_obj_arr);

    wifi_util_dbg_print(WIFI_WEBCONFIG, "%s:%d: Encoding channel scan response object\n", __func__,
        __LINE__);
    if (encode_em_channel_stats_params(neigh_stats, st_obj_arr) != webconfig_error_none) {
        wifi_util_error_print(WIFI_WEBCONFIG, "%s:%d: Failed to encode channel scan response object\n",
            __func__, __LINE__);
        cJSON_Delete(json);
        return webconfig_error_encode;
    }

    str = cJSON_Print(json);

    data->u.encoded.raw = (webconfig_subdoc_encoded_raw_t)calloc(strlen(str) + 1, sizeof(char));
    if (data->u.encoded.raw == NULL) {
        wifi_util_error_print(WIFI_WEBCONFIG, "%s:%d Failed to allocate memory.\n", __func__,
            __LINE__);
        cJSON_free(str);
        cJSON_Delete(json);
        return webconfig_error_encode;
    }

    wifi_util_dbg_print(WIFI_WEBCONFIG, "%s:%d: Encoded data is %s\n", __func__, __LINE__, str);
    memcpy(data->u.encoded.raw, str, strlen(str));
    cJSON_free(str);
    cJSON_Delete(json);

    return webconfig_error_none;
}

webconfig_error_t decode_em_channel_stats_subdoc(webconfig_t *config, webconfig_subdoc_data_t *data)
{
    webconfig_subdoc_t *doc;
    unsigned int i;
    cJSON *json;
    webconfig_subdoc_decoded_data_t *params;
    wifi_util_dbg_print(WIFI_WEBCONFIG, "%s:%d\n", __func__, __LINE__);

    params = &data->u.decoded;
    if (params == NULL) {
        wifi_util_error_print(WIFI_WEBCONFIG, "%s:%d: NULL Pointer\n", __func__, __LINE__);
        return webconfig_error_decode;
    }

    json = data->u.encoded.json;
    if (json == NULL) {
        wifi_util_error_print(WIFI_WEBCONFIG, "%s:%d: NULL json pointer\n", __func__, __LINE__);
        return webconfig_error_decode;
    }

    char *str;
    str = cJSON_Print(json);
    wifi_util_dbg_print(WIFI_WEBCONFIG, "%s:%d: Decoded Str is : %s\n", __func__, __LINE__, str);
    cJSON_free(str);

    doc = &config->subdocs[data->type];

    for (i = 0; i < doc->num_objects; i++) {
        if ((cJSON_GetObjectItem(json, doc->objects[i].name)) == NULL) {
            wifi_util_error_print(WIFI_WEBCONFIG,
                "%s:%d: object:%s not present, validation failed\n", __func__, __LINE__,
                doc->objects[i].name);
            cJSON_Delete(json);
            return webconfig_error_invalid_subdoc;
        }
    }

    channel_scan_response_t **ch_st = (channel_scan_response_t **)&params->collect_stats.stats;
    if (decode_em_channel_stats_object(ch_st, json) != webconfig_error_none) {
        wifi_util_error_print(WIFI_WEBCONFIG, "%s:%d: Failed to decode channel scan response object\n", __func__,
            __LINE__);
        cJSON_Delete(json);
        free(*ch_st);
        return webconfig_error_invalid_subdoc;
    }

    cJSON_Delete(json);

    return webconfig_error_none;
}
#endif
