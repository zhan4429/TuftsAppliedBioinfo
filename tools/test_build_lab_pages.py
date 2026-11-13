"""Regression checks for answer visibility and Markdown tables in lab pages."""
import unittest
from tools.build_lab_pages import render_markdown

class AnswerRenderingTests(unittest.TestCase):
    def test_answer_is_collapsible_and_closed(self):
        _, _, _, body, _ = render_markdown('# Lab\n\nQuestion\n\n<details>\n<summary>Answer</summary>\n\nAnswer text.\n\n</details>\n')
        self.assertIn('<details>', body)
        self.assertNotIn('<details open', body)
        self.assertLess(body.index('<summary>Answer</summary>'), body.index('Answer text.'))
        self.assertLess(body.index('Answer text.'), body.index('</details>'))

    def test_legacy_answer_callout_is_also_collapsed(self):
        _, _, _, body, _ = render_markdown('# Lab\n\n::: Answer\nAnswer text.\n:::\n')
        self.assertIn('<details><summary>Answer</summary>', body)
        self.assertNotIn(':::', body)
        self.assertNotIn('<details open', body)
        self.assertIn('</details>', body)

    def test_unclosed_legacy_callout_stops_build(self):
        with self.assertRaisesRegex(ValueError, 'Unclosed'):
            render_markdown('# Lab\n\n::: Answer\nAnswer text.\n')

    def test_troubleshooting_callout_has_its_own_label(self):
        _, _, _, body, _ = render_markdown('# Lab\n\n::: Troubleshooting\nCheck paths.\n:::\n')
        self.assertIn('<summary>Troubleshooting</summary>', body)

    def test_fence_text_is_allowed_in_a_code_example(self):
        _, _, _, body, _ = render_markdown('# Lab\n\n```text\n::: Answer\n```\n')
        self.assertIn('<pre', body)
        self.assertIn('::: Answer', body)

class TableRenderingTests(unittest.TestCase):
    def test_shell_pipeline_stays_in_one_cell(self):
        _, _, _, body, _ = render_markdown(
            '# Lab\n\n| Task | Command |\n|---|---|\n'
            '| Align | `bwa-mem2 mem reads.fq \\| samtools sort \\| samtools markdup` |\n'
        )
        self.assertEqual(body.count('<th>'), 2)
        self.assertEqual(body.count('<td>'), 2)
        self.assertIn('<code>bwa-mem2 mem reads.fq | samtools sort | samtools markdup</code>', body)

    def test_logical_or_stays_in_one_cell(self):
        _, _, _, body, _ = render_markdown(
            '# Lab\n\n| Task | Command |\n|---|---|\n'
            '| Filter | `QD<2 \\|\\| FS>60` |\n'
        )
        self.assertEqual(body.count('<td>'), 2)
        self.assertIn('<code>QD&lt;2 || FS&gt;60</code>', body)

    def test_escaped_pipe_can_end_a_cell_without_an_outer_separator(self):
        _, _, _, body, _ = render_markdown(
            '# Lab\n\n| Task | Text |\n|---|---|\n'
            '| Symbol | A literal \\|\n'
        )
        self.assertEqual(body.count('<td>'), 2)
        self.assertIn('<td>A literal |</td>', body)


if __name__ == '__main__':
    unittest.main()
