#!/bin/bash

#####
# This script (part 7) this script writes the last (static) part of the required JavaScript and closes the html file
# this script should be called only once.

OUTPUT="$1"
OUTPUT_HTML="$2"

cat << EOF >> ${OUTPUT_HTML}
    }

    function initTabs() {

        var tabs = [];

        // Grab the tab links and content divs from the page
        var liItems = Array.from(document.getElementById('tabs').childNodes).filter(function(item) {
            return item.nodeName === 'LI'
        });

        var tabLinks = liItems.map(function (li) {
            return getFirstChildWithTagName(li, 'A');
        })

        tabLinks.forEach(function (tabLink) {
            tabLink.onclick = showTab;
            tabLink.onfocus = function () {
                this.blur()
            };
            var id = getHash(tabLink.getAttribute('href'));
            var contentDiv = document.getElementById(id);
            tabs.push({
                id: id,
                link: tabLink,
                content: contentDiv
            });
        })

        tabs[0].content.className = 'tabContent';
        tabs[0].link.className = 'selected';
        for (var i = 1; i < tabs.length; i++) {
            tabs[i].link.className = '';
            tabs[i].content.className = 'tabContent hide'
        }

        function showTab() {

            var selectedId = getHash(this.getAttribute('href'));

            tabs.forEach(function (tab) {

                if (tab.id === selectedId) {
                    tab.link.className = 'selected';
                    tab.content.className = 'tabContent';
                } else {
                    tab.link.className = '';
                    tab.content.className = 'tabContent hide';
                }
            });

            igv.visibilityChange();

            // Stop the browser following the link
            return false;
        }


        function getFirstChildWithTagName(element, tagName) {
            for (var i = 0; i < element.childNodes.length; i++) {
                if (element.childNodes[i].nodeName === tagName) return element.childNodes[i];
            }
        }

        function getHash(url) {
            var hashPos = url.lastIndexOf('#');
            return url.substring(hashPos + 1);
        }
    }


</script>


</body>
</html>
EOF

touch "${OUTPUT}"