#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIST="$ROOT/dist"
APP_NAME="scGUI"
DMG_NAME="scGUI"
VERSION="0.2.1"

APP="$DIST/$APP_NAME.app"
if [[ ! -d "$APP" ]]; then
  echo "ERROR: App bundle not found: $APP" >&2
  echo "Run: $ROOT/scripts/make_app.sh" >&2
  exit 2
fi

STAGE="$(mktemp -d "${TMPDIR:-/tmp}/scanwr-dmg.XXXXXX")"
trap 'rm -rf "$STAGE"' EXIT

cp -R "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"

OUT="$DIST/${DMG_NAME}-${VERSION}.dmg"
rm -f "$OUT"

echo "Creating DMG: $OUT"
hdiutil create -volname "$DMG_NAME" -srcfolder "$STAGE" -ov -format UDZO "$OUT" >/dev/null

echo "OK: $OUT"
