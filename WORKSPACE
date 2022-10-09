# ========================================================================================
# File: WORKSPACE
# Brief: Workspace definition file for Bazel
# Author: Rolfe Power
# ========================================================================================
workspace(name = "com_github_rjpower4_symplectic")

# ----------------------------------------------------------------------------------------
# HTTP Retrieval Utility
# ----------------------------------------------------------------------------------------
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

# ----------------------------------------------------------------------------------------
# Google Test
# ----------------------------------------------------------------------------------------
google_test_commit = "609281088cfefc76f9d0ce82e1ff6c30cc3591e5"

google_test_sha256 = "5cf189eb6847b4f8fc603a3ffff3b0771c08eec7dd4bd961bfd45477dd13eb73"

http_archive(
    name = "com_google_googletest",
    sha256 = google_test_sha256,
    strip_prefix = "googletest-{}".format(google_test_commit),
    urls = [
        "https://github.com/google/googletest/archive/{}.zip".format(google_test_commit),
    ],
)