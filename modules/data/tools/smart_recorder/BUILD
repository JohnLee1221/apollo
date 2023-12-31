load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")
load("//tools:cpplint.bzl", "cpplint")

package(default_visibility = ["//visibility:public"])

cc_binary(
    name = "smart_recorder",
    srcs = ["smart_recorder.cc"],
    deps = [
        ":post_record_processor",
        ":realtime_record_processor",
        ":smart_recorder_gflags",
    ],
)

cc_library(
    name = "post_record_processor",
    srcs = ["post_record_processor.cc"],
    hdrs = ["post_record_processor.h"],
    deps = [
        ":channel_pool",
        ":record_processor",
        "//cyber",
    ],
)

cc_library(
    name = "realtime_record_processor",
    srcs = ["realtime_record_processor.cc"],
    hdrs = ["realtime_record_processor.h"],
    deps = [
        ":channel_pool",
        ":record_processor",
        "//cyber",
        "//cyber/tools/cyber_recorder:recorder",
        "//modules/common_msgs/monitor_msgs:smart_recorder_status_cc_proto",
        "//modules/monitor/common:monitor_manager",
    ],
)

cc_library(
    name = "record_processor",
    srcs = ["record_processor.cc"],
    hdrs = ["record_processor.h"],
    deps = [
        ":bumper_crash_trigger",
        ":drive_event_trigger",
        ":emergency_mode_trigger",
        ":hard_brake_trigger",
        ":interval_pool",
        ":regular_interval_trigger",
        ":small_topics_trigger",
        ":swerve_trigger",
        "//cyber",
        "//modules/common/util",
    ],
)

cc_library(
    name = "smart_recorder_gflags",
    srcs = ["smart_recorder_gflags.cc"],
    hdrs = ["smart_recorder_gflags.h"],
    deps = [
        "@com_github_gflags_gflags//:gflags",
    ],
)

cc_library(
    name = "drive_event_trigger",
    srcs = ["drive_event_trigger.cc"],
    hdrs = ["drive_event_trigger.h"],
    deps = [
        ":trigger_base",
        "//modules/common/adapters:adapter_gflags",
        "//modules/common_msgs/basic_msgs:drive_event_cc_proto",
    ],
)

cc_library(
    name = "emergency_mode_trigger",
    srcs = ["emergency_mode_trigger.cc"],
    hdrs = ["emergency_mode_trigger.h"],
    deps = [
        ":trigger_base",
        "//modules/common_msgs/chassis_msgs:chassis_cc_proto",
        "//modules/common/adapters:adapter_gflags",
        "//modules/common_msgs/control_msgs:control_cmd_cc_proto",
    ],
)

cc_library(
    name = "bumper_crash_trigger",
    srcs = ["bumper_crash_trigger.cc"],
    hdrs = ["bumper_crash_trigger.h"],
    deps = [
        ":trigger_base",
        "//modules/common_msgs/chassis_msgs:chassis_cc_proto",
        "//modules/common/adapters:adapter_gflags",
    ],
)

cc_library(
    name = "hard_brake_trigger",
    srcs = ["hard_brake_trigger.cc"],
    hdrs = ["hard_brake_trigger.h"],
    deps = [
        ":trigger_base",
        "//modules/common_msgs/chassis_msgs:chassis_cc_proto",
        "//modules/common/adapters:adapter_gflags",
        "//modules/common_msgs/control_msgs:control_cmd_cc_proto",
    ],
)

cc_library(
    name = "small_topics_trigger",
    srcs = ["small_topics_trigger.cc"],
    hdrs = ["small_topics_trigger.h"],
    deps = [
        ":channel_pool",
        ":trigger_base",
    ],
)

cc_library(
    name = "swerve_trigger",
    srcs = ["swerve_trigger.cc"],
    hdrs = ["swerve_trigger.h"],
    deps = [
        ":trigger_base",
        "//modules/common_msgs/chassis_msgs:chassis_cc_proto",
        "//modules/common/adapters:adapter_gflags",
        "//modules/common_msgs/control_msgs:control_cmd_cc_proto",
    ],
)

cc_library(
    name = "regular_interval_trigger",
    srcs = ["regular_interval_trigger.cc"],
    hdrs = ["regular_interval_trigger.h"],
    deps = [
        ":trigger_base",
    ],
)

cc_library(
    name = "trigger_base",
    srcs = ["trigger_base.cc"],
    hdrs = ["trigger_base.h"],
    deps = [
        ":interval_pool",
        "//cyber",
        "//modules/data/tools/smart_recorder/proto:smart_recorder_triggers_cc_proto",
    ],
)

cc_library(
    name = "interval_pool",
    srcs = ["interval_pool.cc"],
    hdrs = ["interval_pool.h"],
    deps = [
        "//cyber",
        "@com_google_absl//:absl",
    ],
)

cc_library(
    name = "channel_pool",
    srcs = ["channel_pool.cc"],
    hdrs = ["channel_pool.h"],
    deps = [
        "//cyber",
        "//modules/common/adapters:adapter_gflags",
    ],
)

cpplint()
